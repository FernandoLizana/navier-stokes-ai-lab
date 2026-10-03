FROM python:3.11-slim
WORKDIR /workspace
COPY pyproject.toml README.md LICENSE ./
COPY python ./python
RUN pip install --no-cache-dir -e "./python[dev]"
CMD ["pytest", "python/ns_exploration/tests", "-v"]
