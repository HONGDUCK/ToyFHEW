from note_include.elem.Ring  import Ring
from typing import List, Tuple, TypeAlias

RLWEctxt : TypeAlias = Tuple[Ring, Ring]
RLWEpctxt: TypeAlias = List[RLWEctxt]
RGSWctxt : TypeAlias = Tuple[RLWEpctxt, RLWEpctxt]
