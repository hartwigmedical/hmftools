package com.hartwig.hmftools.cider.primer

import com.hartwig.hmftools.cider.VDJSequence
import com.hartwig.hmftools.cider.VJ

data class Primer(val target: String, val name: String, val sequence: String, val vj: VJ)

data class VdjPrimerMatch(val vdj: VDJSequence, val primer: Primer, val index: Int, val numMismatch: Int, val fullVdjSequence: String)