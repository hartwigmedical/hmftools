package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash

interface AnonymizedId {
    val hash: Hash
    val id: Int
}
