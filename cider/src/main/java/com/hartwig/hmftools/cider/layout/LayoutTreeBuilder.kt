package com.hartwig.hmftools.cider.layout

import org.apache.logging.log4j.LogManager
import java.util.*
import kotlin.collections.ArrayList

class LayoutTreeBuilder(inputReads: List<LayoutTree.Read>,
                        minBaseQuality: Int,
                        minOverlapBases: Int,
                        minBaseHighQualCount: Int)
{
    val minBaseQuality: Byte = minBaseQuality.toByte()
    val minOverlapBases: Int = minOverlapBases
    val minBaseHighQualCount: Int = minBaseHighQualCount
    val inputReadList: List<LayoutTree.Read>

    init
    {
        val readDataMutableList = ArrayList<LayoutTree.Read>()
        readDataMutableList.addAll(inputReads)

        // we want to sort them by the end of the sequence in the layout
        // left to right, i.e starting from the root
        readDataMutableList.sortWith(
                Comparator.comparingInt({ r: LayoutTree.Read -> r.layoutPosition + r.sequence.length })
                    .thenComparingDouble({ r: LayoutTree.Read -> r.baseQualities.average() }) // handle the highest quality ones first
                    .thenComparingInt({ r: LayoutTree.Read -> r.layoutPosition })
                    .thenComparing({ r: LayoutTree.Read -> r.readKey.readName }) // just the final catch all
        )

        inputReadList = readDataMutableList
    }

    fun build(): LayoutTree
    {
        sLogger.info("building layout tree from {} reads", inputReadList.size)

        val layoutTree = LayoutTree(minBaseQuality, minOverlapBases)

        val retryReads = ArrayList<LayoutTree.Read>()

        // go through the read data list, and add one by one to the list of clusters
        // if there are multiple clusters that matches, we choose the highest one
        for (read in inputReadList)
        {
            // add it to the layout tree
            if (!layoutTree.tryAddRead(read))
            {
                retryReads.add(read)
            }
        }

        for (read in retryReads)
        {
            layoutTree.tryAddRead(read)
        }

        sLogger.info("layout tree complete, num levels: {}", layoutTree.numLevels)

        return layoutTree
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(LayoutTreeBuilder::class.java)


    }
}