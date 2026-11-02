FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

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
COPY system_diagnostics.py integrated_view.py literature_view.py ./
COPY genomics/ ./genomics/
COPY literature/ ./literature/
COPY wastewater/ ./wastewater/
COPY templates/ ./templates/
COPY examples/ ./examples/

RUN useradd --create-home --shell /usr/sbin/nologin appuser
USER appuser

EXPOSE 10000

# One worker limits memory use on a Free instance; allow small analyses time to finish.
CMD ["sh", "-c", "exec gunicorn --bind 0.0.0.0:${PORT:-10000} --workers 1 --timeout 180 app:app"]
