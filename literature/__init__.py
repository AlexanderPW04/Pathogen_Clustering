"""Optional literature and metadata enrichment (Sprint 4).

The package finds publications linked to uploaded samples, selects short
relevant passages, asks a replaceable extraction provider (Claude, an offline
mock, or a future locally hosted model) for candidate metadata, verifies every
quotation against the retrieved text, and only changes analysis metadata after
a person approves each proposal. The genomic analyses never depend on it.
"""
