package com.hartwig.hmftools.teal.analysers

import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TeloUtils
import java.util.*

// analyse the telo bam and come up with our own score
// of F1, F2, F4
class TelomereReadsAnalyser
{
    private val mFragmentTypeCount: MutableMap<ReadGroup.FragmentType, Int> = EnumMap(ReadGroup.FragmentType::class.java)
    fun getFragmentTypeCount(fragType: ReadGroup.FragmentType): Int
    {
        return mFragmentTypeCount.getOrDefault(fragType, 0)
    }

    fun onReadGroup(readGroup: ReadGroup?)
    {
        // we want to classify this read group
        val fragType = TeloUtils.classifyFragment(readGroup)
        mFragmentTypeCount[fragType] = mFragmentTypeCount.getOrDefault(fragType, 0) + 1
    }
}