package com.hartwig.hmftools.cider.layout

import com.hartwig.hmftools.cider.ReadKey

class TestLayoutRead(
    val readId: String,
    readKey: ReadKey,
    sequence: String,
    baseQualities: ByteArray,
    alignedPosition: Int)
    : ReadLayout.Read(readKey, sequence, baseQualities, alignedPosition)
{
    override fun copy(alignedPosition: Int): ReadLayout.Read
    {
        return TestLayoutRead(readId, readKey, sequence, baseQualities, alignedPosition)
    }
}
