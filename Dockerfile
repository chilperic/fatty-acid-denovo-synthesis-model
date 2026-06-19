FROM python:3.10-slim
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir -e .[dev]
CMD ["bash", "-lc", "pytest -q && python scripts/01_list_variants.py && python scripts/02_load_data.py && python scripts/03_smoke_simulation.py"]
