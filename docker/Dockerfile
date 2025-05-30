# Use Alpine as the base image
FROM alpine:latest

# Install minimal required libraries
RUN apk add --no-cache \
    ca-certificates \
    musl

# Copy the prebuilt static Rust binary
COPY --from=binary /storm /usr/local/bin/storm

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/storm"]