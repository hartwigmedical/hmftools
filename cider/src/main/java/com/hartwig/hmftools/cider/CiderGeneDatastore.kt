package com.hartwig.hmftools.cider

import com.google.common.collect.*
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrConstantDiversityRegion
import org.apache.logging.log4j.LogManager
import java.util.EnumMap
import java.util.stream.Collectors

// we use immutable collections here, data can be accessed by multiple threads
interface ICiderGeneDatastore
{
    fun getAnchorSequenceSet(geneType: VJGeneType): Set<String>
    fun getByAnchorSequence(anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByGeneLocation(genomicLocation: GenomicLocation): ImmutableCollection<VJAnchorTemplate>
    fun getVjAnchorGeneLocations(): ImmutableCollection<VJAnchorGenomeLocation>
    fun getIgConstantDiversityRegions(): ImmutableCollection<IgTcrConstantDiversityRegion>
}

open class CiderGeneDatastore(vjAnchorTemplates: List<VJAnchorTemplate>, igTcrConstantDiversityRegions: List<IgTcrConstantDiversityRegion>) : ICiderGeneDatastore
{
    private val sLogger = LogManager.getLogger(javaClass)

    // all of the data here are immutable, so we access them from multiple threads.
    private val mAnchorSequenceMap: ImmutableMultimap<String, VJAnchorTemplate>
    private val mGeneTypeAnchorSeqMap: ImmutableMap<VJGeneType, ImmutableMultimap<String, VJAnchorTemplate>>
    private val mGeneLocationTemplateMap: ImmutableMultimap<GenomicLocation, VJAnchorTemplate>
    private val mVjAnchorGenomeLocations: ImmutableList<VJAnchorGenomeLocation>
    private val mIgTcrConstantDiversityRegions: ImmutableList<IgTcrConstantDiversityRegion>

    override fun getAnchorSequenceSet(geneType: VJGeneType): ImmutableSet<String>
    {
        val anchorSeqMap = mGeneTypeAnchorSeqMap[geneType]
        return if (anchorSeqMap != null) anchorSeqMap.keySet() else ImmutableSet.of()
    }

    override fun getByAnchorSequence(anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    {
        return mAnchorSequenceMap[anchorSeq]
    }

    override fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    {
        val anchorSeqMap = mGeneTypeAnchorSeqMap[geneType]
        return if (anchorSeqMap != null) anchorSeqMap[anchorSeq] else ImmutableSet.of()
    }

    override fun getByGeneLocation(genomicLocation: GenomicLocation): ImmutableCollection<VJAnchorTemplate>
    {
        return mGeneLocationTemplateMap[genomicLocation]
    }

    override fun getVjAnchorGeneLocations(): ImmutableList<VJAnchorGenomeLocation>
    {
        return mVjAnchorGenomeLocations
    }

    override fun getIgConstantDiversityRegions(): ImmutableCollection<IgTcrConstantDiversityRegion>
    {
        return mIgTcrConstantDiversityRegions
    }

    init
    {
        val anchorSequenceMap: Multimap<String, VJAnchorTemplate> = ArrayListMultimap.create()
        val geneTypeAnchorSeqMap: MutableMap<VJGeneType, Multimap<String, VJAnchorTemplate>> = EnumMap(VJGeneType::class.java)
        val geneLocationVJGeneMap: Multimap<GenomicLocation, VJAnchorTemplate> = ArrayListMultimap.create()
        val vjAnchorGenomeLocationMap: MutableMap<GenomicLocation, VJGeneType> = HashMap()

        // from this we find all the anchor sequence locations and fix them
        for (gene in vjAnchorTemplates)
        {
            if (gene.anchorLocation != null)
            {
                geneLocationVJGeneMap.put(gene.anchorLocation, gene)

                // we want to check that same location cannot be used by more than one VJ type
                val existingVjGeneType: VJGeneType? = vjAnchorGenomeLocationMap[gene.anchorLocation]
                if (existingVjGeneType != null && existingVjGeneType != gene.type)
                {
                    sLogger.error(
                        "gene location: {} is used by multiple gene type: {} and {}",
                        gene.anchorLocation, existingVjGeneType, gene.type
                    )
                    throw RuntimeException("gene location: ${gene.anchorLocation} is used by multiple gene type: ${existingVjGeneType} and ${gene.type}")
                }
                vjAnchorGenomeLocationMap[gene.anchorLocation] = gene.type
            }
            if (gene.anchorSequence.isNotEmpty())
            {
                anchorSequenceMap.put(gene.anchorSequence, gene)
                geneTypeAnchorSeqMap.computeIfAbsent(gene.type) { o: VJGeneType? -> ArrayListMultimap.create() }
                    .put(gene.anchorSequence, gene)
            }
        }
        mAnchorSequenceMap = ImmutableMultimap.copyOf(anchorSequenceMap)

        // copy to immutable, have to convert each entry to immutable version as well
        mGeneTypeAnchorSeqMap = ImmutableMap.copyOf(
            geneTypeAnchorSeqMap.entries.stream().collect(
                Collectors.toMap(
                    { entry -> entry.key },
                    { entry -> ImmutableMultimap.copyOf(entry.value) }
                )))

        mGeneLocationTemplateMap = ImmutableMultimap.copyOf(geneLocationVJGeneMap)

        mVjAnchorGenomeLocations =
            ImmutableList.copyOf(vjAnchorGenomeLocationMap.entries.map({ o -> VJAnchorGenomeLocation(o.value, o.key) }))

        // also the constant region
        mIgTcrConstantDiversityRegions = ImmutableList.copyOf(igTcrConstantDiversityRegions)

        sLogger.info("found {} gene locations", mGeneLocationTemplateMap.keySet().size)
    }
}