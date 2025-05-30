# Start from an empty base image
FROM scratch

# Copy the prebuilt static Rust binary from a prior build stage
COPY --from=binary /storm /usr/local/bin/storm

# Set entrypoint
ENTRYPOINT ["/usr/local/bin/storm"]