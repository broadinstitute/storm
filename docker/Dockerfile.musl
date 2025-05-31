# Build stage
FROM messense/rust-musl-cross:x86_64-musl as builder

# Set the working directory
WORKDIR /usr/src/storm

# Copy the source code
COPY . .

# Add the musl target and build
RUN rustup target add x86_64-unknown-linux-musl && \
    cargo build --release --target x86_64-unknown-linux-musl

# Final stage
FROM scratch

# Copy the prebuilt static binary from the builder stage
COPY --from=builder /usr/src/storm/target/x86_64-unknown-linux-musl/release/storm /storm

# Set entrypoint
ENTRYPOINT ["/storm"] 