FROM python:3.11-slim-bookworm

# One numerical thread on a 0.1-CPU instance, fewer malloc arenas, the largest
# variant-by-sample analysis a 512 MB instance can finish (see genomics/input_inspection.py),
# and the largest uploaded table it reads safely: 1.5 million cells peaked at about 250 MB
# in the comprehensive audit, 4 million at 440-530 MB (genomics/excel_tables.py, file_formats.py).
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    OPENBLAS_NUM_THREADS=1 \
    OMP_NUM_THREADS=1 \
    MALLOC_ARENA_MAX=2 \
    MPLCONFIGDIR=/tmp/matplotlib \
    MAX_OBSERVATION_CELLS=700000 \
    MAX_WORKSHEET_CELLS=1500000

WORKDIR /app

# MAFFT is an external Linux program, not a Python package.
RUN apt-get update \
    && apt-get install -y --no-install-recommends mafft \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app.py run_pipeline.py validate_inputs.py genomic_regions.py ./
COPY environmental_analysis.py ./
COPY geocode_locations.py align_sequences.py compare_genetic_distances.py ./
COPY add_temporal_geographic.py detect_potential_clusters.py ./
COPY system_diagnostics.py integrated_view.py literature_view.py dashboard_guide.py results_layout.py ./
COPY genomics/ ./genomics/
COPY literature/ ./literature/
COPY wastewater/ ./wastewater/
COPY templates/ ./templates/
COPY examples/ ./examples/

RUN useradd --create-home --shell /usr/sbin/nologin appuser
USER appuser

EXPOSE 10000

# One worker limits memory use on a Free instance. Threads keep the site and its health
# check answering while an analysis runs, and a slow analysis is not killed as a stuck worker;
# the application itself allows one analysis at a time.
CMD ["sh", "-c", "exec gunicorn --bind 0.0.0.0:${PORT:-10000} --workers 1 --worker-class gthread --threads 4 --timeout 180 app:app"]
