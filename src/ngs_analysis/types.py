from pandera.typing import Series
import pandera as pa
from typing import Optional


# TODO: joint uniqueness, DNA characters
# currently handled in load_reference_dna
class ReferenceDna(pa.DataFrameModel):
    source: Optional[Series[str]] = pa.Field(nullable=False, coerce=True)
    name: Series[str] = pa.Field(nullable=False)
    reference_dna: Series[str] = pa.Field(nullable=False)


# simulate based on reference sources
# could also be used to generate comparison of expected vs actual results 
class SamplePlan(pa.DataFrameModel):
    sample: Series[str] = pa.Field(nullable=False, coerce=True)
    source: Series[str] = pa.Field(nullable=False, coerce=True)
    coverage: Optional[Series[float]] = pa.Field(nullable=False, coerce=True)


# simulate based on input sequences
class DnaPlan(pa.DataFrameModel):
    sample: Series[str] = pa.Field(nullable=False, coerce=True)
    reference_dna: Series[str] = pa.Field(nullable=False, unique=True)


class Samples(pa.DataFrameModel):
    sample: Series[str] = pa.Field(nullable=False, unique=True, 
                                   coerce=True)
    fastq_name: Series[str] = pa.Field(nullable=False, unique=True, 
                                       coerce=True)


# build this by using config.yaml to parse ReferenceDna
class Designs(pa.DataFrameModel):
    source: Optional[Series[str]] = pa.Field(nullable=False, coerce=True)
    name: Series[str] = pa.Field(nullable=False, coerce=True)


class Candidates(pa.DataFrameModel):
    read_index: Series[int] = pa.Field(nullable=False)
    query: Series[str] = pa.Field(nullable=False)
    reference: Series[str] = pa.Field(nullable=False)
    reference_name: Series[str] = pa.Field(nullable=False)

