"""Storm - A high-performance data processing framework.

This package provides Python bindings for the Storm data processing framework.
Storm is designed for high-performance data processing and analysis tasks.
"""

import storm.storm as storm  # noqa: F403


def version() -> str:
    """Get the current version of the Storm package.

    Returns:
        str: The version string of the currently installed Storm package.
    """
    return storm._version()
