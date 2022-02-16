package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator
import kotlin.math.abs

// combine tumor and germline data

class SampleBreakEndFinder(config: BreakEndParams)
{
    private val mConfig: BreakEndParams = config
    private val mBreakEndEvidence = ArrayList<TelomericBreakEndEvidence>()

    fun run()
    {
        // we first collect all the read groups

        // find the tumour breakends
        val tumorBreakEndFinder = TelbamBreakEndFinder(mConfig, java.io.File(mConfig.tumorTelbamFile))
        tumorBreakEndFinder.run()

        // for the germline we also include the break ends we found in the tumor together to investigate
        val germlineBreakEndFinder = TelbamBreakEndFinder(mConfig, java.io.File(mConfig.germlineTelbamFile),
            tumorBreakEndFinder.breakEnds.map({ telomericBreakEnd -> telomericBreakEnd.key }))
        germlineBreakEndFinder.run()

        matchTumorGermlineBreakEnds(tumorBreakEndFinder.breakEnds.toMutableList(),
            germlineBreakEndFinder.breakEnds.toMutableList(), mConfig.markDuplicateDistance)

        markDuplicateBreakEnds()
        BreakEndTsvWriter.writeBreakEnds(mConfig.outputFile, mBreakEndEvidence)
    }

    fun matchTumorGermlineBreakEnds(tumorBreakEnds: List<TelomericBreakEndSupport>,
                                    germlineBreakEnds: List<TelomericBreakEndSupport>,
                                    breakpointMergeDistance: Int)
    {
        mBreakEndEvidence.clear()

        // we sort them both first
        val comparator = Comparator.comparing(TelomericBreakEndSupport::type)
            .thenComparing(TelomericBreakEndSupport::chromosome)
            .thenComparingInt(TelomericBreakEndSupport::position)

        // first sort all the breakends
        // we want to sort it such that position is sorted last, and then we can use that
        // to compare if positions are close by, we merge them
        val tumorItr = tumorBreakEnds.toMutableList().apply({ sortWith(comparator) }).listIterator()
        val germlineItr = germlineBreakEnds.toMutableList().apply({ sortWith(comparator) }).listIterator()

        while (tumorItr.hasNext() && germlineItr.hasNext())
        {
            val tumorBreakEnd =  tumorItr.next()
            val germlineBreakEnd = germlineItr.next()

            var evidence: TelomericBreakEndEvidence?

            // see if they are close to each other, and maybe combine them
            if (tumorBreakEnd.type == germlineBreakEnd.type &&
                tumorBreakEnd.chromosome == germlineBreakEnd.chromosome &&
                abs(tumorBreakEnd.position - germlineBreakEnd.position) <= breakpointMergeDistance)
            {
                // same position, we can combine them
                evidence = TelomericBreakEndEvidence(tumorSupport = tumorBreakEnd, germlineSupport = germlineBreakEnd)
            }
            else if (comparator.compare(tumorBreakEnd, germlineBreakEnd) > 0)
            {
                // germline one is smaller in the sorting order
                evidence = TelomericBreakEndEvidence(germlineSupport = germlineBreakEnd)
                tumorItr.previous()
            }
            else
            {
                // tumor one is smaller in the sorting order
                evidence = TelomericBreakEndEvidence(tumorSupport = tumorBreakEnd)
                germlineItr.previous()
            }

            mBreakEndEvidence.add(evidence)
        }

        // now add any remaining ones, only one of the lists should have anything left
        germlineItr.forEachRemaining({ o -> mBreakEndEvidence.add(TelomericBreakEndEvidence(germlineSupport = o))})
        tumorItr.forEachRemaining({ o-> mBreakEndEvidence.add(TelomericBreakEndEvidence(tumorSupport = o))})

        // sort it by chromosome position
        mBreakEndEvidence.sortWith(
            Comparator.comparing({ o: TelomericBreakEndEvidence -> o.chromosome() }, ContigComparator.INSTANCE)
            .thenComparingInt({ o -> o.position() }))
    }

    fun markDuplicateBreakEnds()
    {
        data class DupCandidate
            (val breakEnd: TelomericBreakEndSupport,
            val inTumor: Boolean)

        fun handleDuplicates(duplicateBreakEnds: List<DupCandidate>)
        {
            if (duplicateBreakEnds.size <= 1)
            {
                return
            }
            // find the main breakend, and mark the rest as duplicates
            // we prefer the main breakend to be from tumor, only if it is not found
            // in tumor we then try to use a germline one
            var mainBreakEndCandidates = duplicateBreakEnds.filter({o -> o.inTumor})

            if (mainBreakEndCandidates.isEmpty())
                mainBreakEndCandidates = duplicateBreakEnds

            val mainBreakEnd = mainBreakEndCandidates.maxByOrNull({ o -> o.breakEnd.numSplitReads() })

            for (candidate in duplicateBreakEnds)
            {
                if (candidate != mainBreakEnd)
                {
                    candidate.breakEnd.key.isDuplicate = true
                }
            }
        }

        // now we go through this list, and for ones that are close to each other, we mark them
        // as duplicate
        var currentBreakEnd: TelomericBreakEnd? = null
        val duplicateBreakEnds = ArrayList<DupCandidate>()

        for (evidence in mBreakEndEvidence)
        {
            val mainBreakEnd = evidence.mainBreakEnd()

            if (mainBreakEnd == null)
                continue

            if (currentBreakEnd != null && currentBreakEnd.type === mainBreakEnd.type &&
                currentBreakEnd.chromosome == evidence.chromosome() &&
                abs(currentBreakEnd.position - evidence.position()) <= mConfig.markDuplicateDistance
            )
            {
                duplicateBreakEnds.add(DupCandidate(breakEnd = mainBreakEnd, inTumor = evidence.tumorSupport != null))
            }
            else
            {
                // this is a new breakend, so finish work on the previous one
                handleDuplicates(duplicateBreakEnds)
                duplicateBreakEnds.clear()
                duplicateBreakEnds.add(DupCandidate(breakEnd = mainBreakEnd, inTumor = evidence.tumorSupport != null))

                currentBreakEnd = mainBreakEnd.key
            }
        }
        handleDuplicates(duplicateBreakEnds)
    }
}

