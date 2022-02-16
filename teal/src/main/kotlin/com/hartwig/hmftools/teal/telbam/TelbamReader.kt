package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.common.genome.chromosome.ContigComparator
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.teal.ReadGroup
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager

// read the telbam file and gives the read groups
class TelbamReader(
    telbamFile: java.io.File,
    excludedRegions: List<GenomeRegion> = emptyList(),
    includedRegions: List<GenomeRegion>? = null)
{
    private val mLogger = LogManager.getLogger(TelbamReader::class.java)
    private val mTelbamFile = telbamFile
    private val mExcludedGenomeRegions = excludedRegions
    private val mIncludedGenomeRegions = includedRegions

    // the distance which we deem the same breakpoint
    private val mReadGroups: MutableMap<String, ReadGroup> = HashMap()
    val readGroups: Map<String, ReadGroup> get() { return mReadGroups }

    // create into read groups
    fun read()
    {
        mLogger.info("processing telbam: {}", mTelbamFile)
        val factory = SamReaderFactory.makeDefault()
        val samReader = factory.open(mTelbamFile)
        samReader.iterator().use({ iterator ->
            while (iterator.hasNext())
            {
                val record = iterator.next()
                processReadRecord(record)
            }
        })
    }

    private fun processReadRecord(record: SAMRecord)
    {
        var readGroup: ReadGroup? = mReadGroups[record.readName]
        if (readGroup == null)
        {
            // cache if new
            readGroup = ReadGroup(record.readName)
            mReadGroups[record.readName] = readGroup
        }
        readGroup.acceptRead(record)

        if (readGroup.isComplete())
        {
            // now we see if all the reads are in the excluded region, if so just remove it
            if (shouldExcludeReadGroup(readGroup))
            {
                mReadGroups.remove(readGroup.name)
            }
        }

        assert(readGroup.invariant())
    }

    fun shouldExcludeReadGroup(readGroup: ReadGroup): Boolean
    {
        if (mIncludedGenomeRegions != null)
        {
            // if none is included then we exclude read group
            return readGroup.allReads.none({ r ->

                var isIncluded = false
                // check if this read should be included
                for (includedRegion in mIncludedGenomeRegions)
                {
                    if (r.alignmentStart <= includedRegion.end() &&
                        r.alignmentEnd >= includedRegion.start() &&
                        ContigComparator.INSTANCE.compare(r.referenceName, includedRegion.chromosome()) == 0)
                    {
                        isIncluded = true
                        break
                    }
                }
                // if not explicitly included then it is excluded
                isIncluded
            })
        }
        else
        {
            // if all reads are excluded then we exclude read group
            return readGroup.allReads.all({ r ->
                var isExcluded = false
                for (excludedRegion in mExcludedGenomeRegions)
                {
                    if (r.alignmentStart <= excludedRegion.end() &&
                        r.alignmentEnd >= excludedRegion.start() &&
                        ContigComparator.INSTANCE.compare(r.referenceName, excludedRegion.chromosome()) == 0
                    )
                    {
                        isExcluded = true
                        break
                    }
                }
                isExcluded
            })
        }
    }
}
