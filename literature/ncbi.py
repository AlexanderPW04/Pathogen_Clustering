"""NCBI E-utilities: follow record links from sample identifiers to PubMed.

Only NCBI's own record links are used (for example nuccore_pubmed), so every
association found here is traceable to an NCBI record. Documentation:
https://www.ncbi.nlm.nih.gov/books/NBK25497/
"""

from literature.http_client import ApiClient

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "bmi461_genomic_surveillance"
SEARCH_DATABASES = {"Nucleotide": ("nuccore", "{}[Accession]"), "BioSample": ("biosample", "{}"),
                    "BioProject": ("bioproject", "{}"), "SRA": ("sra", "{}")}


class NcbiClient:
    def __init__(self, api_key="", email="", max_requests=80, session=None):
        # Without an API key NCBI allows 3 requests per second; with one, 10.
        interval = 0.11 if api_key else 0.34
        self.client = ApiClient("NCBI E-utilities", [EUTILS], interval, max_requests, session=session)
        self.base = {"tool": TOOL_NAME, "retmode": "json"}
        if api_key:
            self.base["api_key"] = api_key
        if email:
            self.base["email"] = email

    def call(self, name, **params):
        return self.client.get(EUTILS + name, {**self.base, **params})

    def uids(self, database, term):
        result = self.call("esearch.fcgi", db=database, term=term, retmax=5)
        return list(result.get("esearchresult", {}).get("idlist", []))

    def links(self, database, uid, target):
        """[(linkname, [ids])] from one record to a target database."""

        result = self.call("elink.fcgi", dbfrom=database, db=target, id=uid)
        found = []
        for linkset in result.get("linksets", []):
            for group in linkset.get("linksetdbs", []) or []:
                found.append((group.get("linkname", ""), [str(item) for item in group.get("links", [])]))
        return found

    def pubmed_for_record(self, kind, value):
        """Publications linked from one sample record: [(pmid, relationship, evidence)]."""

        database, template = SEARCH_DATABASES[kind]
        results = []
        for uid in self.uids(database, template.format(value))[:1]:
            direct = "study" if kind == "BioProject" else "direct"
            for linkname, pmids in self.links(database, uid, "pubmed"):
                for pmid in pmids[:5]:
                    results.append((pmid, direct, f"NCBI {database} record {value} links PubMed {pmid} ({linkname})."))
            if kind != "BioProject":
                for _, projects in self.links(database, uid, "bioproject"):
                    for project in projects[:2]:
                        for linkname, pmids in self.links("bioproject", project, "pubmed"):
                            for pmid in pmids[:5]:
                                results.append((pmid, "study",
                                                f"NCBI {database} record {value} belongs to BioProject uid "
                                                f"{project}, which links PubMed {pmid} ({linkname})."))
        return results

    def pmid_for(self, kind, value):
        """Convert a PMCID or DOI to a PMID with a PubMed search."""

        field = {"PMCID": "pmcid", "DOI": "doi"}[kind]
        found = self.uids("pubmed", f'"{value}"[{field}]')
        return found[0] if len(found) == 1 else None

    def summaries(self, pmids):
        """Title, date, and identifiers for PubMed records."""

        if not pmids:
            return {}
        result = self.call("esummary.fcgi", db="pubmed", id=",".join(pmids)).get("result", {})
        summaries = {}
        for pmid in pmids:
            record = result.get(str(pmid))
            if not isinstance(record, dict) or "error" in record:
                continue
            ids = {item.get("idtype"): item.get("value") for item in record.get("articleids", [])}
            summaries[str(pmid)] = {
                "Title": record.get("title", ""),
                "Publication_Date": record.get("pubdate", ""),
                "PMID": str(pmid), "PMCID": ids.get("pmc", ""), "DOI": ids.get("doi", ""),
                "Source_URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        return summaries
