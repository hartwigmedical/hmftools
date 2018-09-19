package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash

data class HashId(override val hash: Hash, override val id: Int) : AnonymizedId {
    companion object {
        operator fun invoke(hashString: String, id: Int) = HashId(Hash(hashString), id)
    }
}
