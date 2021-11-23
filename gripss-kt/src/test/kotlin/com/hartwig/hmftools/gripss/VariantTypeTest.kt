package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.gripsskt.*
import org.junit.Assert.assertEquals
import org.junit.Test

class VariantTypeTest {

    @Test
    fun testMobileElementInsertion(){
        isNotMobileElementInsert("T", "TTTTTTTTTT.")
        isNotMobileElementInsert("C", "CTTTTTTTTT.")
        isMobileElementInsertion("C", "CTTTTTTTTTT.")
        isNotMobileElementInsert("C", "CTTTTTTTTTG.")

        isNotMobileElementInsert("C", "CTTTTTTTTTGTTTTTT.")
        isMobileElementInsertion("C", "CTTTTTTTTTGTTTTTTT.")

        isNotMobileElementInsert("A", ".AAAAAAAAAA")
        isNotMobileElementInsert("C", ".AAAAAAAAAC")
        isMobileElementInsertion("C", ".AAAAAAAAAAC")
        isNotMobileElementInsert("C", ".GAAAAAAAAAC")

        isNotMobileElementInsert("C", ".AAAAAAAGAAGAAAAAC")
        isMobileElementInsertion("C", ".AAAAAAAAAAGAAAAAAC")
    }

    private fun isMobileElementInsertion(ref: String, alt: String) {
        isMobileElementInsert(true, ref, alt)

        if (alt.endsWith(".")) {
            isMobileElementInsert(true, ref, alt.replace(".", "[1:1000["))

            val reverseAlt = "." + alt.replace(".","'")
            isMobileElementInsert(false, ref, reverseAlt)
            isMobileElementInsert(false, ref, reverseAlt.replace(".", "]1:1000]"))
        }

        if (alt.startsWith(".'")) {
            isMobileElementInsert(true, ref, alt.replace(".", "]1:1000]"))

            val reverseAlt = alt.replace(".","'") + "."
            isMobileElementInsert(false, ref, reverseAlt)
            isMobileElementInsert(false, ref, reverseAlt.replace(".", "[1:1000["))
        }
    }

    private fun isNotMobileElementInsert(ref: String, alt: String) {
        isMobileElementInsert(false, ref, alt)
    }

    private fun isMobileElementInsert(expected: Boolean, ref: String, alt: String) {
        assertEquals(expected, VariantType.create("1", 100, ref, alt).isMobileElementInsertion())
    }

    @Test
    fun testExtraColonInContig() {
        val start = VariantType.create("1", 100, "A", "AA[HLA:10:100:1000[")
        assertEquals(Translocation("A", "A","HLA:10:100", 1000, 1, -1), start)
    }

    @Test
    fun testSingle() {
        val sgl = VariantType.create("1", 100, "A", "TAAA.");
        assertEquals(Single("T", "AAA", 1), sgl)
    }

    @Test
    fun testBnd() {
        val start = VariantType.create("1", 100, "A", "AAT[2:110[")
        assertEquals(Translocation("A", "AT","2", 110, 1, -1), start)

        val end = VariantType.create("2", 110, "T", "]1:100]AAT")
        assertEquals(Translocation("T", "AA","1", 100, -1, 1), end)
    }

    @Test
    fun testDup() {
        val start = VariantType.create("1", 110, "A", "AAT[1:100[")
        assertEquals(Duplication("A", "AT","1", 100, 1, -1, 10), start)

        val end = VariantType.create("1", 100, "T", "]1:110]AAT")
        assertEquals(Duplication("T", "AA","1", 110, -1, 1, 10), end)
    }

    @Test
    fun testDel() {
        val start = VariantType.create("1", 100, "A", "AAT[1:110[")
        assertEquals(Deletion("A", "AT","1", 110, 1, -1, 10), start)

        val end = VariantType.create("1", 110, "T", "]1:100]AAT")
        assertEquals(Deletion("T", "AA","1", 100, -1, 1, 10), end)
    }

    @Test
    fun testIns() {
        val start = VariantType.create("1", 100, "A", "AAT[1:101[")
        assertEquals(Insertion("A", "AT","1", 101, 1, -1, 1), start)

        val end = VariantType.create("1", 101, "T", "]1:100]AAT")
        assertEquals(Insertion("T", "AA","1", 100, -1, 1, 1), end)
    }

    @Test
    fun testToString() {
        assertToString(".GATAC")
        assertToString("GATAC.")
        assertToString("CAT[1:3123[")
        assertToString("CAT]1:3123]")
        assertToString("[1:3123[CAT")
        assertToString("]1:3123]CAT")
    }

    private fun assertToString(expected: String) {
        val variant = VariantType.create("1", 100, "A", expected)
        assertEquals(expected, variant.toString())
    }
}