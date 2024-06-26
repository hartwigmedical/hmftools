package com.hartwig.hmftools.teal.telbam

import com.hartwig.hmftools.common.bam.SamRecordUtils
import com.hartwig.hmftools.teal.TealUtils
import htsjdk.samtools.SamReader
import org.apache.logging.log4j.LogManager

class PartitionReader(private val bamPartition: BamPartition,
                      private val writer: BamRecordWriter,
                      private val incompleteReadNames: MutableSet<String>)
{
    fun readPartition(samReader: SamReader)
    {
        bamPartition.iterator(samReader).use { itr ->
            var numReads = 0L
            for (read in itr)
            {
                if (read.hasAttribute(SamRecordUtils.CONSENSUS_READ_ATTRIBUTE))
                {
                    // we want to ignore consensus read in TEAL
                    continue
                }

                val hasTeloContent = TealUtils.hasTelomericContent(read.readString)
                if (hasTeloContent || incompleteReadNames.contains(read.readName))
                {
                    writer.processReadRecord(read, hasTeloContent)
                    ++numReads
                }
            }
            logger.info("processed partition({}) num reads({})", bamPartition, numReads)
        }
    }

    companion object
    {
        private val logger = LogManager.getLogger(javaClass)
    }
}