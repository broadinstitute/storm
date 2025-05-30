# Use a multi-architecture base image
FROM --platform=$BUILDPLATFORM alpine:latest

# Install minimal required libraries
RUN apk add --no-cache \
    ca-certificates \
    musl

# Copy the appropriate binary based on architecture
ARG TARGETARCH
COPY --from=binary /storm /usr/local/bin/storm

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/storm"]