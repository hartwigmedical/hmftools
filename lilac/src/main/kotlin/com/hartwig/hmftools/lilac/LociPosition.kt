package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion
import com.hartwig.hmftools.common.genome.region.Strand

class LociPosition(private val transcripts: List<HmfTranscriptRegion>) {

    fun position(codingLoci: Int): List<Int> {
        val result = mutableListOf<Int>()
        for (transcript in transcripts) {
            if (transcript.strand() == Strand.FORWARD) {
                result.add(forwardPosition(codingLoci, transcript))
            } else {
                result.add(reversePosition(codingLoci, transcript))
            }
        }
        return result
    }

    fun reversePosition(codingLoci: Int, transcript: HmfTranscriptRegion): Int {
        require(transcript.strand() == Strand.REVERSE)

        val codingStart = transcript.codingStart();
        val codingEnd = transcript.codingEnd();

        var currentLoci = 0
        for (exon in transcript.exome().reversed()) {
            if (exon.end() >= codingStart && exon.start() <= codingEnd) {
                val exonStartPosition = Math.max(codingStart, exon.start())
                val exonEndPosition = Math.min(codingEnd, exon.end())
                val exonLength = (exonEndPosition - exonStartPosition + 1).toInt()

                val exonStartLoci = currentLoci
                val exonEndLoci = exonStartLoci + exonLength - 1

                if (codingLoci in exonStartLoci..exonEndLoci) {
                    return (exonEndPosition - codingLoci + exonStartLoci).toInt()
                }

                currentLoci = exonEndLoci + 1
            }
        }

        return -1
    }

    fun forwardPosition(codingLoci: Int, transcript: HmfTranscriptRegion): Int {
        require(transcript.strand() == Strand.FORWARD)

        val codingStart = transcript.codingStart();
        val codingEnd = transcript.codingEnd();


        var currentLoci = 0
        for (exon in transcript.exome()) {
            if (exon.end() >= codingStart && exon.start() <= codingEnd) {
                val exonStartPosition = Math.max(codingStart, exon.start())
                val exonEndPosition = Math.min(codingEnd, exon.end())
                val exonLength = (exonEndPosition - exonStartPosition + 1).toInt()

                val exonStartLoci = currentLoci
                val exonEndLoci = exonStartLoci + exonLength - 1

                if (codingLoci in exonStartLoci..exonEndLoci) {
                    return (exonStartPosition + codingLoci - exonStartLoci).toInt()
                }

                currentLoci = exonEndLoci + 1
            }
        }

        return -1
    }

}