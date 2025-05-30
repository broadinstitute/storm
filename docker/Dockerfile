# Use a minimal base image with required libraries
FROM debian:bookworm-slim

# Install minimal required libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    libc6 \
    && rm -rf /var/lib/apt/lists/*

# Copy the prebuilt static Rust binary
COPY --from=binary /storm /usr/local/bin/storm

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/storm"]