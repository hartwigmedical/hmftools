package com.hartwig.hmftools.cider;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.eclipse.collections.api.collection.ImmutableCollection;
import org.eclipse.collections.api.factory.Maps;
import org.eclipse.collections.api.factory.Sets;
import org.eclipse.collections.api.list.ImmutableList;
import org.eclipse.collections.api.map.ImmutableMap;
import org.eclipse.collections.api.map.MutableMap;
import org.eclipse.collections.api.multimap.ImmutableMultimap;
import org.eclipse.collections.api.multimap.MutableMultimap;
import org.eclipse.collections.api.set.SetIterable;
import org.eclipse.collections.impl.list.mutable.FastList;
import org.eclipse.collections.impl.map.mutable.UnifiedMap;
import org.eclipse.collections.impl.multimap.list.FastListMultimap;
import org.jetbrains.annotations.NotNull;

public class CiderGeneDataLoader implements CiderGeneDatastore
{
    private static final Logger sLogger = LogManager.getLogger(CiderGeneDataLoader.class);

    // all of the data here are immutable, so we access them from multiple threads.
    private final ImmutableMultimap<String, VJAnchorTemplate> mAnchorSequenceMap;

    private final ImmutableMap<VJGeneType, ImmutableMultimap<String, VJAnchorTemplate>> mGeneTypeAnchorSeqMap;

    private final ImmutableMultimap<VJAnchorReferenceLocation, VJAnchorTemplate> mGeneLocationVJGeneMap;

    private final ImmutableList<IgTcrConstantRegion> mIgTcrConstantRegions;

    @NotNull
    @Override
    public SetIterable<String> getAnchorSequenceSet(@NotNull VJGeneType geneType)
    {
        ImmutableMultimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.keySet() : Sets.immutable.empty();
    }

    @NotNull
    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorSequence(@NotNull String anchorSeq)
    {
        return mAnchorSequenceMap.get(anchorSeq);
    }

    @NotNull
    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorSequence(@NotNull VJGeneType geneType, @NotNull String anchorSeq)
    {
        ImmutableMultimap<String, VJAnchorTemplate> anchorSeqMap = mGeneTypeAnchorSeqMap.get(geneType);
        return anchorSeqMap != null ? anchorSeqMap.get(anchorSeq) : Sets.immutable.empty();
    }

    @NotNull
    @Override
    public ImmutableCollection<VJAnchorTemplate> getByAnchorGeneLocation(@NotNull VJAnchorReferenceLocation vjAnchorReferenceLocation)
    {
        return mGeneLocationVJGeneMap.get(vjAnchorReferenceLocation);
    }

    @NotNull
    @Override
    public SetIterable<VJAnchorReferenceLocation> getVJAnchorReferenceLocations()
    {
        return mGeneLocationVJGeneMap.keySet();
    }

    @NotNull
    @Override
    public ImmutableCollection<IgTcrConstantRegion> getIgConstantRegions() { return mIgTcrConstantRegions; }

    public CiderGeneDataLoader(RefGenomeVersion refGenomeVersion, String ensemblDataDir) throws IOException
    {
        List<VJAnchorTemplate> vjAnchorTemplates = loadAnchorTemplateTsv(refGenomeVersion);

        MutableMultimap<String, VJAnchorTemplate> anchorSequenceMap = new FastListMultimap<>();
        MutableMap<VJGeneType, MutableMultimap<String, VJAnchorTemplate>> geneTypeAnchorSeqMap = new UnifiedMap<>();
        MutableMultimap<VJAnchorReferenceLocation, VJAnchorTemplate> geneLocationVJGeneMap = new FastListMultimap<>();

        // from this we find all the anchor sequence locations and fix them
        for (VJAnchorTemplate gene : vjAnchorTemplates)
        {
            if (gene.getAnchorLocation() != null)
            {
                geneLocationVJGeneMap.put(new VJAnchorReferenceLocation(gene.getType().getVj(), gene.getAnchorLocation()), gene);
            }

            if (!gene.getAnchorSequence().isEmpty())
            {
                anchorSequenceMap.put(gene.getAnchorSequence(), gene);
                geneTypeAnchorSeqMap.computeIfAbsent(gene.getType(), o -> new FastListMultimap<>())
                    .put(gene.getAnchorSequence(), gene);
            }
        }

        mIgTcrConstantRegions = loadConstantRegionGenes(refGenomeVersion, ensemblDataDir);

        mAnchorSequenceMap = anchorSequenceMap.toImmutable();

        // copy to immutable, have to convert each entry to immutable version as well
        mGeneTypeAnchorSeqMap = Maps.immutable.ofMap(geneTypeAnchorSeqMap.entrySet().stream().collect(
                        Collectors.toMap(Map.Entry::getKey, e -> e.getValue().toImmutable())));
        mGeneLocationVJGeneMap = geneLocationVJGeneMap.toImmutable();

        sLogger.info("found {} gene locations", mGeneLocationVJGeneMap.keySet().size());
    }

    private static List<VJAnchorTemplate> loadAnchorTemplateTsv(RefGenomeVersion refGenomeVersion) throws IOException
    {
        String resourcePath;

        if (refGenomeVersion.is37())
        {
            resourcePath = "igtcr_anchor.37.tsv";
        }
        else if (refGenomeVersion.is38())
        {
            resourcePath = "igtcr_anchor.38.tsv";
        }
        else
        {
            throw new IllegalArgumentException("unknown ref genome version: " + refGenomeVersion);
        }
        List<VJAnchorTemplate> VJAnchorTemplateList = new ArrayList<>();

        java.io.InputStream tsvStream = CiderGeneDataLoader.class.getClassLoader().getResourceAsStream(resourcePath);
        if (tsvStream == null)
        {
            sLogger.error("unable to find resource file: {}", resourcePath);
            throw new RuntimeException("unable to find resource file: " + resourcePath);
        }

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(tsvStream)))
        {
            CSVFormat format = CSVFormat.Builder.create()
                .setDelimiter('\t')
                .setRecordSeparator('\n')
                .setHeader().setSkipHeaderRecord(true) // use first line header as column names
                .build();
            Iterable<CSVRecord> records = format.parse(reader);
    
            for (CSVRecord record : records)
            {
                String id = record.get("id");
                String name = record.get("gene");
                String allele = record.get("allele");
                int posStart = Integer.parseInt(record.get("posStart"));
                int posEnd = Integer.parseInt(record.get("posEnd"));
                String anchorSequence = record.get("anchorSequence");
                String chromosome = record.get("chr");

                String strandStr = record.get("strand");
                Strand strand = null;
                if (strandStr.equals("+"))
                    strand = Strand.FORWARD;
                else if (strandStr.equals("-"))
                    strand = Strand.REVERSE;
                int anchorStart = Integer.parseInt(record.get("anchorStart"));
                int anchorEnd = Integer.parseInt(record.get("anchorEnd"));
                String sequence = record.get("sequence");

                GeneLocation geneLocation = null;
                GeneLocation anchorLocation = null;

                if (!chromosome.isEmpty())
                {
                    chromosome = refGenomeVersion.versionedChromosome(chromosome);

                    if (posStart <= 0 || posEnd <= 0)
                    {
                        throw new RuntimeException("chromosome exist but pos start or pos end invalid");
                    }

                    if (strand == null)
                    {
                        throw new RuntimeException("chromosome exist but strand invalid");
                    }

                    geneLocation = new GeneLocation(chromosome, posStart, posEnd, strand);

                    if (anchorStart >= 0 && anchorEnd >= 0)
                        anchorLocation = new GeneLocation(chromosome, anchorStart, anchorEnd, strand);
                }

                VJAnchorTemplate VJAnchorTemplate = new VJAnchorTemplate(
                        id, name, allele, geneLocation, sequence, anchorSequence, anchorLocation);
                VJAnchorTemplateList.add(VJAnchorTemplate);
            }
        }
        return VJAnchorTemplateList;
    }

    private ImmutableList<IgTcrConstantRegion> loadConstantRegionGenes(RefGenomeVersion refGenomeVersion, String ensemblDataDir)
    {
        final EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDataDir, refGenomeVersion);
        boolean ensemblLoadOk = ensemblDataCache.load(true);

        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed");
            throw new RuntimeException("Ensembl data cache load failed");
        }

        FastList<IgTcrConstantRegion> igTcrConstantRegions = new FastList<>();

        // find all the constant region genes
        Map<String,List<GeneData>> chrGeneDataMap = ensemblDataCache.getChrGeneDataMap();
        List<GeneData> geneDataList = chrGeneDataMap.values().stream().flatMap(Collection::stream).collect(Collectors.toList());

        // find every gene that is constant region
        for (GeneData geneData : geneDataList)
        {
            if (geneData.GeneName.length() <= 4)
                continue;

            String geneNamePrefix = geneData.GeneName.substring(0, 4);
            IgTcrConstantRegion.Type igConstantRegionType;

            try
            {
                igConstantRegionType = IgTcrConstantRegion.Type.valueOf(geneNamePrefix);
            }
            catch (IllegalArgumentException ignored)
            {
                continue;
            }

            GeneLocation geneLocation = new GeneLocation(geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd,
                    geneData.forwardStrand() ? Strand.FORWARD : Strand.REVERSE);

            igTcrConstantRegions.add(new IgTcrConstantRegion(igConstantRegionType, geneLocation));

            sLogger.info("found constant region gene: {}, type: {}, location: {}",
                    geneData.GeneName, igConstantRegionType, geneLocation);
        }

        return igTcrConstantRegions.toImmutable();
    }
}
