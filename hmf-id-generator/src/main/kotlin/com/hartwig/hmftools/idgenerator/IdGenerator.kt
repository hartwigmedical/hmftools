package com.hartwig.hmftools.idgenerator

import org.bouncycastle.jcajce.provider.digest.SHA3
import org.bouncycastle.util.encoders.Hex

class IdGenerator(private val password: String) {

    fun hash(plaintext: String): String {
        val sha3 = SHA3.Digest256()
        sha3.update((plaintext + password).toByteArray())
        return Hex.toHexString(sha3.digest())
    }
}
