package com.hartwig.hmftools.lilackt

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import junit.framework.Assert.assertEquals
import junit.framework.Assert.assertTrue
import org.junit.Test

class LociPositionTest {

    private val transcripts = HmfGenePanelSupplier.allGenesMap37()
    private val aTranscript = transcripts[LilacApplication.HLA_A]!!
    private val bTranscript = transcripts[LilacApplication.HLA_B]!!
    private val hlaTranscripts = listOf(aTranscript, bTranscript, transcripts[LilacApplication.HLA_C]!!)
    private val victim = LociPosition(hlaTranscripts)

    @Test
    fun testA() {
        val minLoci = victim.nucelotideLoci(aTranscript.codingStart().toInt())
        assertTrue(minLoci >= 0)

        val maxLoci = victim.nucelotideLoci(aTranscript.codingEnd().toInt())

        for (i in minLoci..maxLoci) {
            val position = victim.position(i, aTranscript)
            val loci = victim.nucelotideLoci(position)
            assertEquals(i, loci)
        }
    }

    @Test
    fun testB() {
        val minLoci = victim.nucelotideLoci(bTranscript.codingEnd().toInt())
        val maxLoci = victim.nucelotideLoci(bTranscript.codingStart().toInt())
        assertTrue(minLoci >= 0)
        assertTrue(maxLoci > minLoci)

        for (i in minLoci..maxLoci) {
            val position = victim.position(i, bTranscript)
            val loci = victim.nucelotideLoci(position)
            assertEquals(i, loci)
        }
    }

}