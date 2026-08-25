"""Common predictor contract."""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

@dataclass(frozen=True)
class PredictionResult:
    model: str
    prediction: str
    probability: float | None = None
    extra: dict[str, Any] = field(default_factory=dict)

class PredictorUnavailable(RuntimeError):
    pass

class BasePredictor(ABC):
    name = "Predictor"
    @abstractmethod
    def predict(self, sequence: str) -> PredictionResult: ...
    @abstractmethod
    def availability(self) -> tuple[bool, str]: ...
