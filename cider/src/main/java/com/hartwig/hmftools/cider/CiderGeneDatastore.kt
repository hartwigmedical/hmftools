package com.hartwig.hmftools.cider

import com.google.common.collect.*
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrConstantDiversityRegion
import com.hartwig.hmftools.cider.genes.VJAnchorGenomeLocation
import com.hartwig.hmftools.cider.genes.VJAnchorTemplate
import com.hartwig.hmftools.cider.genes.VJGeneType
import org.apache.logging.log4j.LogManager
import java.util.EnumMap

// we use immutable collections here, data can be accessed by multiple threads
interface ICiderGeneDatastore
{
    fun getAnchorSequenceSet(geneType: VJGeneType): List<String>
    fun getByAnchorSequence(anchorSeq: String): List<VJAnchorTemplate>
    fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): List<VJAnchorTemplate>
    fun getByGeneLocation(genomicLocation: GenomicLocation): List<VJAnchorTemplate>
    fun getVjAnchorGeneLocations(): List<VJAnchorGenomeLocation>
    fun getIgConstantDiversityRegions(): List<IgTcrConstantDiversityRegion>
}

open class CiderGeneDatastore(vjAnchorTemplates: List<VJAnchorTemplate>, igTcrConstantDiversityRegions: List<IgTcrConstantDiversityRegion>) : ICiderGeneDatastore
{
    private val sLogger = LogManager.getLogger(javaClass)

    // all of the data here are immutable, so we access them from multiple threads.
    private val mAnchorSequenceMap: Map<String, List<VJAnchorTemplate>>
    private val mGeneTypeAnchorSeqMap: Map<VJGeneType, Map<String, List<VJAnchorTemplate>>>
    private val mGeneTypeAnchorSequences: Map<VJGeneType, List<String>>
    private val mGeneLocationTemplateMap: Map<GenomicLocation, List<VJAnchorTemplate>>
    private val mVjAnchorGenomeLocations: List<VJAnchorGenomeLocation>
    private val mIgTcrConstantDiversityRegions: List<IgTcrConstantDiversityRegion>
    
    override fun getAnchorSequenceSet(geneType: VJGeneType): List<String>
    {
        return mGeneTypeAnchorSequences[geneType] ?: emptyList()
    }

    override fun getByAnchorSequence(anchorSeq: String): List<VJAnchorTemplate>
    {
        return mAnchorSequenceMap[anchorSeq] ?: emptyList()
    }

    override fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): List<VJAnchorTemplate>
    {
        val anchorSeqMap = mGeneTypeAnchorSeqMap[geneType]
        return anchorSeqMap?.get(anchorSeq) ?: emptyList()
    }

    override fun getByGeneLocation(genomicLocation: GenomicLocation): List<VJAnchorTemplate>
    {
        return mGeneLocationTemplateMap[genomicLocation] ?: emptyList()
    }

    override fun getVjAnchorGeneLocations(): List<VJAnchorGenomeLocation>
    {
        return mVjAnchorGenomeLocations
    }

    override fun getIgConstantDiversityRegions(): List<IgTcrConstantDiversityRegion>
    {
        return mIgTcrConstantDiversityRegions
    }

    init
    {
        // Care is taken here to compute these fields in a deterministic manner to allow reproducible results between runs.
        // Note that hash-based structures do not give consistent ordering.

        val anchorSequenceMap: MutableMap<String, MutableList<VJAnchorTemplate>> = HashMap()
        val geneTypeAnchorSeqMap: MutableMap<VJGeneType, MutableMap<String, MutableList<VJAnchorTemplate>>> = EnumMap(VJGeneType::class.java)
        val geneLocationVJGeneMap: MutableMap<GenomicLocation, MutableList<VJAnchorTemplate>> = HashMap()
        val vjAnchorGenomeLocationMap: MutableMap<GenomicLocation, VJGeneType> = HashMap()
        val vjAnchorGenomeLocations = ArrayList<VJAnchorGenomeLocation>()

        // from this we find all the anchor sequence locations and fix them
        for (gene in vjAnchorTemplates)
        {
            if (gene.anchorLocation != null)
            {
                geneLocationVJGeneMap.computeIfAbsent(gene.anchorLocation) { ArrayList() }.add(gene)

                // we want to check that same location cannot be used by more than one VJ type
                val existingVjGeneType: VJGeneType? = vjAnchorGenomeLocationMap[gene.anchorLocation]
                if (existingVjGeneType == null)
                {
                    vjAnchorGenomeLocationMap[gene.anchorLocation] = gene.type
                    vjAnchorGenomeLocations.add(VJAnchorGenomeLocation(gene.type, gene.anchorLocation))
                }
                else if (existingVjGeneType != gene.type)
                {
                    sLogger.error(
                        "gene location: {} is used by multiple gene type: {} and {}",
                        gene.anchorLocation, existingVjGeneType, gene.type
                    )
                    throw RuntimeException("gene location: ${gene.anchorLocation} is used by multiple gene type: ${existingVjGeneType} and ${gene.type}")
                }
            }
            if (gene.anchorSequence.isNotEmpty())
            {
                anchorSequenceMap.computeIfAbsent(gene.anchorSequence) { ArrayList() }.add(gene)
                geneTypeAnchorSeqMap.computeIfAbsent(gene.type) { HashMap() }
                    .computeIfAbsent(gene.anchorSequence) { ArrayList() }.add(gene)
            }
        }

        mAnchorSequenceMap = ImmutableMap.copyOf(anchorSequenceMap.mapValues
            {entry -> ImmutableList.copyOf(entry.value)})

        // copy to immutable, have to convert each entry to immutable version as well
        mGeneTypeAnchorSeqMap = ImmutableMap.copyOf(
            geneTypeAnchorSeqMap.mapValues
                    { entry -> ImmutableMap.copyOf(entry.value.mapValues
                            { o -> ImmutableList.copyOf(o.value) })
                    })
        // Precompute this because requesting the keys of mGeneTypeAnchorSeqMap will give inconsistent ordering, which affects downstream results.
        mGeneTypeAnchorSequences = ImmutableMap.copyOf(
            mGeneTypeAnchorSeqMap.mapValues { entry -> ImmutableList.copyOf(entry.value.keys.sorted()) })

        mGeneLocationTemplateMap = ImmutableMap.copyOf(geneLocationVJGeneMap.mapValues
            {entry -> ImmutableList.copyOf(entry.value)})

        mVjAnchorGenomeLocations = ImmutableList.copyOf(vjAnchorGenomeLocations)

        // also the constant region
        mIgTcrConstantDiversityRegions = ImmutableList.copyOf(igTcrConstantDiversityRegions)

        sLogger.info("found {} gene locations", mGeneLocationTemplateMap.keys.size)
    }
}