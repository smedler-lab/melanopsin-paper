FROM mcr.microsoft.com/devcontainers/python:1-3.10-bookworm

# Install R (r-base pulls in a recent R release from the Debian bookworm repos).
# Drop the base image's yarn apt source first: its signing key is expired/unavailable
# and it isn't needed for this project, but its presence breaks `apt-get update`.
RUN rm -f /etc/apt/sources.list.d/yarn.list

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
