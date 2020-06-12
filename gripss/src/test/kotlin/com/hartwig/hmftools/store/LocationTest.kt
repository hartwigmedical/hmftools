package com.hartwig.hmftools.store

import com.hartwig.hmftools.bedpe.Breakend
import org.junit.Assert
import org.junit.Test

class LocationTest {

    @Test
    fun testBEALN() {
        val line = "X:52729547|+|608M|0,X:52786405|-|608M|,X:48209476|+|608M|".split(",")
        val breakends = line.map { Breakend.fromBealn(it) }

        Assert.assertEquals(Breakend("X", 52729547, 52729547, 1), breakends[0])
        Assert.assertEquals(Breakend("X", 52786405, 52786405, -1), breakends[1])
        Assert.assertEquals(Breakend("X", 48209476, 48209476, 1), breakends[2])
    }
}
