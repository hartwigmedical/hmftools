package com.hartwig.hmftools.bedpe

import com.hartwig.hmftools.gripss.store.LocationString
import org.junit.Assert.assertEquals
import org.junit.Test

class LocationStringTest {

    @Test
    fun testFromString() {
        assertEquals(LocationString("chr1", 10000), LocationString("chr1:10000"))
        assertEquals(LocationString("HLA-A*01:04N", 1737), LocationString("HLA-A*01:04N:1737"))
        assertEquals(LocationString("HLA:10:100", 1000), LocationString("HLA:10:100:1000"))
        assertEquals(LocationString("kraken:taxid_67082_NC_032111.1", 8387), LocationString("kraken:taxid_67082_NC_032111.1:8387"))
    }

}