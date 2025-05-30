# Use a multi-architecture base image
FROM --platform=$BUILDPLATFORM alpine:latest

# Install minimal required libraries
RUN apk add --no-cache \
    ca-certificates \
    musl

# Copy the prebuilt static Rust binary
COPY --from=binary /storm /usr/local/bin/storm

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/storm"]