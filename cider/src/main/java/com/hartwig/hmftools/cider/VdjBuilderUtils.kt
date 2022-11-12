package com.hartwig.hmftools.cider

import org.apache.logging.log4j.LogManager

object VdjBuilderUtils
{
    private val sLogger = LogManager.getLogger(javaClass)

    data class VjLayoutOverlap(
        val vLayoutOffset: Int, // offset in overlap
        val jLayoutOffset: Int, // offset in overlap
        val overlap: Int,
        val highQualMatchBases: Int)

    // we try to overlap them. Note that it is possible for the V and J layout to overlap in any
    // direction
    //
    // Normally we have?
    // ----------------++++++++++++++++++++++         v layout
    //                 ++++++++++++++++++++++------   j layout
    //                 |--- overlap (>0) ---|
    //
    // it is also possible that
    //                 ++++++++++++++++++++++------   v layout
    // ----------------++++++++++++++++++++++         j layout
    //                 |--- overlap (<0) ---|
    fun findVjOverlap(vLayoutSeq: String, jLayoutSeq: String, minOverlappedBases: Int) : VjLayoutOverlap?
    {
        var highQualMatchBases: Int = 0

        // first we try sliding the J seq along the V seq
        for (i in 0 until vLayoutSeq.length - minOverlappedBases)
        {
            var seqMatch = true
            highQualMatchBases = 0
            val overlapSize = Math.min(vLayoutSeq.length - i, jLayoutSeq.length)

            // check for overlap
            for (j in 0 until overlapSize)
            {
                val vIndex: Int = i + j
                val jIndex: Int = j

                val vBase = vLayoutSeq[vIndex]
                val jBase = jLayoutSeq[jIndex]

                val bothHighQual = (vBase != 'N' && jBase != 'N')

                if (bothHighQual)
                {
                    if (vBase == jBase)
                    {
                        ++highQualMatchBases
                    }
                    else
                    {
                        seqMatch = false
                        break
                    }
                }
            }

            if (seqMatch && highQualMatchBases > minOverlappedBases)
            {
                // found overlap
                return VjLayoutOverlap(0, i, overlapSize, highQualMatchBases)
            }
        }

        // second we try sliding the V seq along the J seq
        for (i in 0 until jLayoutSeq.length - minOverlappedBases)
        {
            var seqMatch = true
            highQualMatchBases = 0
            val overlapSize = Math.min(jLayoutSeq.length - i, vLayoutSeq.length)

            // check for overlap
            for (j in 0 until overlapSize)
            {
                val vIndex: Int = j
                val jIndex: Int = i + j

                val vBase = vLayoutSeq[vIndex]
                val jBase = jLayoutSeq[jIndex]

                val bothHighQual = (vBase != 'N' && jBase != 'N')

                if (bothHighQual)
                {
                    if (vBase == jBase)
                    {
                        ++highQualMatchBases
                    }
                    else
                    {
                        seqMatch = false
                        break
                    }
                }
            }

            if (seqMatch && highQualMatchBases > minOverlappedBases)
            {
                // found overlap
                return VjLayoutOverlap(i, 0, overlapSize, highQualMatchBases)
            }
        }

        return null
    }

    // we already worked out where the overlaps are, we now want to work out the
    // distance between V and J
    fun calcVJDistance(vSeq: String, jSeq: String, vPosInVSeq: Int, jPosInJSeq: Int, overlap: Int) : Int
    {
        // v layout:   TGC-GAATACC-CACATCCTGA-G
        // j layout:            CC-CACATCCTGA-GAGTG-TCAGA
        //                 |_____|            |___|
        //            V anchor(align)    J anchor(align)

        // v layout:    GC-GAATACC-CACATCCTGA-GAGTG-TCAGA
        // j layout: GATGC-GAATACC-CACATCCTGA-GAGTG-T
        //                 |_____|            |___|
        //            V anchor(align)    J anchor(align)
        //                        |----------|

        return 0
    }
}