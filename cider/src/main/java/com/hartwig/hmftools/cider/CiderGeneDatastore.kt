package com.hartwig.hmftools.cider

import org.apache.logging.log4j.LogManager
import org.eclipse.collections.api.collection.ImmutableCollection
import org.eclipse.collections.api.factory.Lists
import org.eclipse.collections.api.factory.Maps
import org.eclipse.collections.api.factory.Sets
import org.eclipse.collections.api.list.ImmutableList
import org.eclipse.collections.api.map.ImmutableMap
import org.eclipse.collections.api.map.MutableMap
import org.eclipse.collections.api.multimap.ImmutableMultimap
import org.eclipse.collections.api.multimap.MutableMultimap
import org.eclipse.collections.api.set.SetIterable
import org.eclipse.collections.impl.map.mutable.UnifiedMap
import org.eclipse.collections.impl.multimap.list.FastListMultimap
import java.util.stream.Collectors

// we use immutable collections here, data can be accessed by multiple threads
interface ICiderGeneDatastore
{
    fun getAnchorSequenceSet(geneType: VJGeneType): SetIterable<String>
    fun getByAnchorSequence(anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByGeneLocation(genomeRegionStrand: GenomeRegionStrand): ImmutableCollection<VJAnchorTemplate>
    fun getVjAnchorGeneLocations(): ImmutableCollection<VJAnchorGenomeLocation>
    fun getIgConstantRegions(): ImmutableCollection<IgTcrConstantRegion>
}

open class CiderGeneDatastore(vjAnchorTemplates: List<VJAnchorTemplate>, igTcrConstantRegions: List<IgTcrConstantRegion>) : ICiderGeneDatastore
{
    private val sLogger = LogManager.getLogger(javaClass)

    // all of the data here are immutable, so we access them from multiple threads.
    private val mAnchorSequenceMap: ImmutableMultimap<String, VJAnchorTemplate>
    private val mGeneTypeAnchorSeqMap: ImmutableMap<VJGeneType, ImmutableMultimap<String, VJAnchorTemplate>>
    private val mGeneLocationTemplateMap: ImmutableMultimap<GenomeRegionStrand, VJAnchorTemplate>
    private val mVjAnchorGenomeLocations: ImmutableList<VJAnchorGenomeLocation>
    private val mIgTcrConstantRegions: ImmutableList<IgTcrConstantRegion>

    override fun getAnchorSequenceSet(geneType: VJGeneType): SetIterable<String>
    {
        val anchorSeqMap = mGeneTypeAnchorSeqMap[geneType]
        return if (anchorSeqMap != null) anchorSeqMap.keySet() else Sets.immutable.empty()
    }

    override fun getByAnchorSequence(anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    {
        return mAnchorSequenceMap[anchorSeq]
    }

    override fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    {
        val anchorSeqMap = mGeneTypeAnchorSeqMap[geneType]
        return if (anchorSeqMap != null) anchorSeqMap[anchorSeq] else Sets.immutable.empty()
    }

    override fun getByGeneLocation(genomeRegionStrand: GenomeRegionStrand): ImmutableCollection<VJAnchorTemplate>
    {
        return mGeneLocationTemplateMap[genomeRegionStrand]
    }

    override fun getVjAnchorGeneLocations(): ImmutableList<VJAnchorGenomeLocation>
    {
        return mVjAnchorGenomeLocations
    }

    override fun getIgConstantRegions(): ImmutableCollection<IgTcrConstantRegion>
    {
        return mIgTcrConstantRegions
    }

    init
    {
        val anchorSequenceMap: MutableMultimap<String, VJAnchorTemplate> = FastListMultimap()
        val geneTypeAnchorSeqMap: MutableMap<VJGeneType, MutableMultimap<String, VJAnchorTemplate>> = UnifiedMap()
        val geneLocationVJGeneMap: MutableMultimap<GenomeRegionStrand, VJAnchorTemplate> = FastListMultimap()
        val vjAnchorGenomeLocationMap: MutableMap<GenomeRegionStrand, VJGeneType> = UnifiedMap()

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
                geneTypeAnchorSeqMap.computeIfAbsent(gene.type) { o: VJGeneType? -> FastListMultimap() }
                    .put(gene.anchorSequence, gene)
            }
        }
        mAnchorSequenceMap = anchorSequenceMap.toImmutable()

        // copy to immutable, have to convert each entry to immutable version as well
        mGeneTypeAnchorSeqMap = Maps.immutable.ofMap(
            geneTypeAnchorSeqMap.entries.stream().collect(
                Collectors.toMap(
                    { entry -> entry.key },
                    { entry -> entry.value.toImmutable() }
                )))

        mGeneLocationTemplateMap = geneLocationVJGeneMap.toImmutable()

        mVjAnchorGenomeLocations =
            Lists.immutable.fromStream(vjAnchorGenomeLocationMap.entries.stream().map({ o -> VJAnchorGenomeLocation(o.value, o.key) }))

        // also the constant region
        mIgTcrConstantRegions = Lists.immutable.ofAll(igTcrConstantRegions)

        sLogger.info("found {} gene locations", mGeneLocationTemplateMap.keySet().size())
    }
}