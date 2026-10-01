# pleioFDR container image (linux/amd64). Build:
#   docker build --platform linux/amd64 -t pleiofdr .
# The image needs no network access at run time; data are bind-mounted, see README "Containers".

# python:3.13-slim-bookworm (multi-arch index), pinned for reproducible builds
ARG BASE=python:3.13-slim-bookworm@sha256:2325bb286ec344af3e5898cc224b5844e2707ac6e26b1632516fd3edc84a5e26

FROM ghcr.io/astral-sh/uv:0.12.18 AS uv

FROM ${BASE} AS builder
COPY --from=uv /uv /usr/local/bin/uv
ENV UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_DOWNLOADS=never \
    UV_PYTHON=/usr/local/bin/python3.13 \
    UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy
WORKDIR /build
# locked third-party dependencies first, so code changes do not reinstall them
COPY pyproject.toml uv.lock ./
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-install-project
COPY README.md LICENSE ./
COPY src ./src
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev --no-editable

FROM ${BASE}
ARG VERSION=2.0.0
LABEL org.opencontainers.image.title="pleiofdr" \
      org.opencontainers.image.description="Pleiotropy-informed conditional and conjunctional false discovery rate" \
      org.opencontainers.image.source="https://github.com/precimed/pleiofdr" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.version="${VERSION}"
COPY --from=builder /opt/venv /opt/venv
COPY config_default.txt config_template.txt /opt/pleiofdr/share/
# numba and matplotlib cache folders are chosen at run time by pleiofdr (pleiofdr/_caches.py),
# so the image works for any user id and with a read-only root filesystem
ENV PATH=/opt/venv/bin:$PATH \
    MPLBACKEND=Agg \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1
RUN python -c "import pleiofdr.analysis, pleiofdr.fuma.combine, pleiofdr.fuma.novelty" \
    && rm -rf /tmp/* /root/.cache /root/.config \
    && useradd --create-home --uid 1000 pleiofdr
USER pleiofdr
WORKDIR /data
ENTRYPOINT ["pleiofdr"]
CMD ["--help"]
