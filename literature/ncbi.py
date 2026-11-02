"""NCBI E-utilities: follow record links from sample identifiers to PubMed.

Only NCBI's own record links are used (for example nuccore_pubmed), so every
association found here is traceable to an NCBI record. Documentation:
https://www.ncbi.nlm.nih.gov/books/NBK25497/

Requests are batched: ESummary accepts a comma-separated list of accession
versions, and ELink with repeated `id` parameters returns one link set per
record, so ten accessions need a handful of requests, not ten of each.
Shared BioProjects and BioSamples are looked up once.
"""

from literature.http_client import ApiClient

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "bmi461_genomic_surveillance"
SEARCH_DATABASES = {"Nucleotide": ("nuccore", "{}[Accession]"), "BioSample": ("biosample", "{}"),
                    "BioProject": ("bioproject", "{}"), "SRA": ("sra", "{}")}
BATCH = 100  # identifiers per ESummary / ELink request


def chunks(values, size=BATCH):
    values = list(values)
    return [values[start:start + size] for start in range(0, len(values), size)]


class NcbiClient:
    def __init__(self, api_key="", email="", max_requests=80, session=None):
        # Without an API key NCBI allows 3 requests per second; with one, 10.
        interval = 0.11 if api_key else 0.34
        self.client = ApiClient("NCBI E-utilities", [EUTILS], interval, max_requests, session=session)
        self.base = {"tool": TOOL_NAME, "retmode": "json"}
        if api_key:
            self.base["api_key"] = api_key  # sent to NCBI only; never logged or shown
        if email:
            self.base["email"] = email

    def call(self, name, **params):
        return self.client.get(EUTILS + name, {**self.base, **params})

    def uids(self, database, term):
        result = self.call("esearch.fcgi", db=database, term=term, retmax=5)
        return list(result.get("esearchresult", {}).get("idlist", []))

    def links(self, database, uid, target):
        """[(linkname, [ids])] from one record to a target database."""

        return self.links_by_id(database, [uid], target).get(str(uid), [])

    def links_by_id(self, database, uids, target):
        """{uid: [(linkname, [ids])]} for many records, one link set per record."""

        found = {}
        for group in chunks(dict.fromkeys(str(uid) for uid in uids)):
            result = self.call("elink.fcgi", dbfrom=database, db=target, id=group)
            for linkset in result.get("linksets", []):
                ids = [str(item) for item in linkset.get("ids", [])]
                pairs = [(item.get("linkname", ""), [str(link) for link in item.get("links", [])])
                         for item in linkset.get("linksetdbs", []) or []]
                for uid in ids or group[:1]:
                    found.setdefault(uid, []).extend(pairs)
        return found

    def nucleotide_uids(self, accessions):
        """{accession as given: uid} from ESummary, which accepts accession versions in batches."""

        wanted = {str(value).split(".")[0].upper(): str(value) for value in accessions}
        found = {}
        for group in chunks(dict.fromkeys(str(value) for value in accessions)):
            result = self.call("esummary.fcgi", db="nuccore", id=",".join(group)).get("result", {})
            for uid in result.get("uids", []):
                record = result.get(uid, {})
                for name in [record.get("accessionversion", ""), record.get("caption", "")]:
                    key = str(name).split(".")[0].upper()
                    if key in wanted:
                        found[wanted[key]] = str(uid)
        return found

    def pubmed_for_records(self, kind, values):
        """{value: [(pmid, relationship, evidence)]} for many identifiers of one kind.

        Nucleotide records: direct PubMed links, links through each record's
        BioSample (one specimen, so sample-specific), and links through its
        BioProject (a study, so study-level). Other kinds are looked up one
        value at a time (they are usually shared by many samples).
        """

        database, template = SEARCH_DATABASES[kind]
        if kind == "Nucleotide":
            uid_of = self.nucleotide_uids(values)
        else:
            uid_of = {}
            for value in dict.fromkeys(values):
                uids = self.uids(database, template.format(value))[:1]
                if uids:
                    uid_of[value] = uids[0]
        uids = list(dict.fromkeys(uid_of.values()))
        results = {value: [] for value in values}
        if not uids:
            return results
        relationship = "study" if kind == "BioProject" else "direct"
        direct = self.links_by_id(database, uids, "pubmed")
        biosamples = self.links_by_id(database, uids, "biosample") if kind == "Nucleotide" else {}
        projects = self.links_by_id(database, uids, "bioproject") if kind != "BioProject" else {}
        sample_ids = list(dict.fromkeys(link for pairs in biosamples.values() for _, links in pairs
                                        for link in links[:2]))
        project_ids = list(dict.fromkeys(link for pairs in projects.values() for _, links in pairs
                                         for link in links[:2]))
        sample_pubmed = self.links_by_id("biosample", sample_ids, "pubmed") if sample_ids else {}
        project_pubmed = self.links_by_id("bioproject", project_ids, "pubmed") if project_ids else {}
        for value, uid in uid_of.items():
            found = results.setdefault(value, [])
            for linkname, pmids in direct.get(uid, []):
                found += [(pmid, relationship, f"NCBI {database} record {value} links PubMed {pmid} ({linkname}).")
                          for pmid in pmids[:5]]
            for _, links in biosamples.get(uid, []):
                for biosample in links[:2]:
                    for linkname, pmids in sample_pubmed.get(biosample, []):
                        found += [(pmid, "direct", f"NCBI {database} record {value} has BioSample uid {biosample}, "
                                                   f"which links PubMed {pmid} ({linkname}).") for pmid in pmids[:5]]
            for _, links in projects.get(uid, []):
                for project in links[:2]:
                    for linkname, pmids in project_pubmed.get(project, []):
                        found += [(pmid, "study", f"NCBI {database} record {value} belongs to BioProject uid "
                                                  f"{project}, which links PubMed {pmid} ({linkname}).")
                                  for pmid in pmids[:5]]
            results[value] = list(dict.fromkeys(found))
        return results

    def pubmed_for_record(self, kind, value):
        """Publications linked from one sample record: [(pmid, relationship, evidence)]."""

        return self.pubmed_for_records(kind, [value]).get(value, [])

    def biosample_accessions_for(self, nucleotides):
        """{nucleotide accession: [BioSample accessions]}: other names for the same samples.

        The ELink requests repeat ones made by pubmed_for_records and are
        answered from the cache; one ESummary names all BioSamples.
        """

        uid_of = self.nucleotide_uids(nucleotides)
        links = self.links_by_id("nuccore", list(uid_of.values()), "biosample") if uid_of else {}
        sample_of = {value: [link for _, ids in links.get(uid, []) for link in ids[:2]]
                     for value, uid in uid_of.items()}
        unique = list(dict.fromkeys(link for ids in sample_of.values() for link in ids))
        names = {}
        for group in chunks(unique):
            result = self.call("esummary.fcgi", db="biosample", id=",".join(group)).get("result", {})
            names.update({uid: str(result[uid]["accession"]) for uid in group
                          if isinstance(result.get(uid), dict) and result[uid].get("accession")})
        return {value: [names[uid] for uid in dict.fromkeys(ids) if uid in names]
                for value, ids in sample_of.items()}

    def biosample_accessions(self, nucleotide):
        return self.biosample_accessions_for([nucleotide]).get(nucleotide, [])

    def pmid_for(self, kind, value):
        """Convert a PMCID or DOI to a PMID with a PubMed search."""

        field = {"PMCID": "pmcid", "DOI": "doi"}[kind]
        found = self.uids("pubmed", f'"{value}"[{field}]')
        return found[0] if len(found) == 1 else None

    def summaries(self, pmids):
        """Title, dates, and identifiers for PubMed records (one request per 100 PMIDs)."""

        summaries = {}
        for group in chunks(dict.fromkeys(str(pmid) for pmid in pmids)):
            result = self.call("esummary.fcgi", db="pubmed", id=",".join(group)).get("result", {})
            for pmid in group:
                record = result.get(str(pmid))
                if not isinstance(record, dict) or "error" in record:
                    continue
                ids = {item.get("idtype"): item.get("value") for item in record.get("articleids", [])}
                summaries[str(pmid)] = {
                    "Title": record.get("title", ""),
                    "Publication_Date": record.get("pubdate", ""),
                    "Epub_Date": record.get("epubdate", ""),
                    "PMID": str(pmid), "PMCID": ids.get("pmc", ""), "DOI": ids.get("doi", ""),
                    "Source_URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                }
        return summaries
