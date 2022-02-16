package com.hartwig.hmftools.teal.breakend

import com.hartwig.hmftools.teal.telbam.TelbamReader
import org.apache.logging.log4j.LogManager

class TelbamBreakEndFinder(config: BreakEndParams, telbamFile: java.io.File, extraCandidateBreakEnds: Collection<TelomericBreakEnd> = ArrayList())
{
    private val mLogger = LogManager.getLogger(TelbamBreakEndFinder::class.java)
    private val mConfig = config
    private val mTelbamFile = telbamFile

    // extra breakends that we want to investigate in addition to the ones
    // found by split read analyser
    private val mExtraCandidateBreakEnds = extraCandidateBreakEnds

    var breakEnds: List<TelomericBreakEndSupport> = ArrayList()

    fun run()
    {
        mLogger.info("reading telbam file: {}", mTelbamFile.path)

        // first we want to find the read groups
        val telbamReader = TelbamReader(mTelbamFile, mConfig.excludedGenomeRegions, mConfig.includedGenomeRegions)
        telbamReader.read()

        val readGroups = telbamReader.readGroups.values

        // now we got the read groups
        val splitReadAnalyser = CandidateBreakEndFinder(
            mConfig,
            extraCandidateBreakEnds = mExtraCandidateBreakEnds)

        splitReadAnalyser.processReadGroups(readGroups)

        // now we got all the potential break ends, we send them off to the
        val breakEndSupportCounter = BreakEndSupportCounter(mConfig.refGenomeVersion, mConfig.telomereMatchThreshold)
        breakEnds = breakEndSupportCounter.countBreakEndSupports(splitReadAnalyser.potentialBreakEnds, readGroups)

        // we perform a clean up to remove break ends that have 0 support
        breakEnds = breakEnds.filter({ b -> b.fragments.isNotEmpty() }).toList()
        // breakEnds = breakEnds.filter({ b -> b.longestSplitReadAlignLength >= mConfig.MinSplitReadAlignLength }).toList()
    }
}


