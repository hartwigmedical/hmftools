package htsjdk.samtools

import junit.framework.Assert.assertEquals
import org.junit.Test

class SAMUtilsExtraTest {
    @Test
    fun testSoftClipsBeforeAndAfterMatch() {
        val cigar =  TextCigarCodec.decode("10S130M10S")

        val result = SAMUtilsExtra.getAlignmentBlocksWithSoftClips(cigar, 100, "TypeName")
        assertEquals(3, result.size)
        assertAlignment(result[0], 90, 1, 10)
        assertAlignment(result[1], 100, 11, 130)
        assertAlignment(result[2], 230, 141, 10)
    }

    @Test
    fun testSoftClipsInterruptedByIndelBeforeAndAfterMatch() {
        val cigar = TextCigarCodec.decode("10S1D130M1D10S")

        val result = SAMUtilsExtra.getAlignmentBlocksWithSoftClips(cigar, 100, "TypeName")
        assertEquals(1, result.size)
        assertAlignment(result[0], 101, 11, 130)
    }

    private fun assertAlignment(victim: AlignmentBlock, expectedReferenceStart: Int, expectedReadStart: Int, expectedLength: Int) {
        assertEquals(expectedReferenceStart, victim.referenceStart)
        assertEquals(expectedReadStart, victim.readStart)
        assertEquals(expectedLength, victim.length)
    }

}