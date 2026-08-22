from abc import ABC, abstractmethod
from typing import Any

from deconvatac.data import DeconvolutionInput, DeconvolutionResult


class BaseDeconvolver(ABC):
    """Base class for method adapters using the shared input/result contracts."""

    method_name: str

    def __init__(self, **kwargs: Any):
        self.config = kwargs

    @abstractmethod
    def run(self, data: DeconvolutionInput) -> DeconvolutionResult:
        """Run deconvolution and return standardized results."""
