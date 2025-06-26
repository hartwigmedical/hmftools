package com.hartwig.hmftools.linx.visualiser;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.linx.visualiser.SvVisualiser.VIS_LOGGER;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_COPY_NUMBER_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_FUSIONS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_GENE_EXONS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_PROTEIN_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_LINKS_FILE;
import static com.hartwig.hmftools.linx.visualiser.file.VisDataWriter.COHORT_VIS_SVS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.linx.visualiser.data.VisCopyNumbers;
import com.hartwig.hmftools.linx.visualiser.data.VisExons;
import com.hartwig.hmftools.linx.visualiser.data.VisProteinDomains;
import com.hartwig.hmftools.linx.visualiser.data.VisSegments;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

public class SampleData
{
    public final List<VisSegment> Segments;
    public final List<VisSvData> SvData;
    public final List<VisCopyNumber> CopyNumbers;
    public final List<VisProteinDomain> ProteinDomains;
    public final List<VisFusion> Fusions;
    public final List<VisGeneExon> Exons;

    public final List<AmberBAF> AmberBAFs = Lists.newArrayList();
    public final List<CobaltRatio> CobaltRatios = Lists.newArrayList();
    public final List<PurpleSegment> PurpleSegments = Lists.newArrayList();

    private final VisualiserConfig mConfig;

    public SampleData(final VisualiserConfig config) throws Exception
    {
        mConfig = config;

        boolean isGermline = config.IsGermline;

        if(!isGermline
        && !Files.exists(Paths.get(VisSvData.generateFilename(mConfig.SampleDataDir, mConfig.Sample, false)))
        && Files.exists(Paths.get(VisSvData.generateFilename(mConfig.SampleDataDir, mConfig.Sample, true))))
        {
            isGermline = true;
        }

        final String svDataFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_SVS_FILE : VisSvData.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        final String linksFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_LINKS_FILE : VisSegment.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        final String cnaFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_COPY_NUMBER_FILE : VisCopyNumber.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        final String geneExonFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_GENE_EXONS_FILE : VisGeneExon.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        final String proteinFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_PROTEIN_FILE : VisProteinDomain.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        final String fusionFile = mConfig.UseCohortFiles ?
                mConfig.SampleDataDir + COHORT_VIS_FUSIONS_FILE : VisFusion.generateFilename(mConfig.SampleDataDir, mConfig.Sample, isGermline);

        List<VisSvData> svData = VisSvData.read(svDataFile).stream().filter(x -> matchOnSampleId(x.SampleId)).collect(toList());

        boolean svDataFiltered = false;

        if(!mConfig.SpecificRegions.isEmpty())
        {
            List<VisSvData> svsInRegions = svData.stream()
                    .filter(x -> mConfig.SpecificRegions.stream()
                            .anyMatch(y -> y.containsPosition(x.ChrStart, x.PosStart) || y.containsPosition(x.ChrEnd, x.PosEnd)))
                    .collect(toList());

            if(mConfig.ClusterIds.isEmpty())
            {
                // if clusters have not been specified, then add these to the set to be plotted
                svsInRegions.forEach(x -> addClusterId(x.ClusterId));
            }
            else
            {
                // otherwise only show SVs for the specified clusters which are also in the specified regions
                svData = svsInRegions;
            }

            svDataFiltered = true;
        }

        if(!mConfig.ClusterIds.isEmpty())
        {
            svData = svData.stream().filter(x -> mConfig.ClusterIds.contains(x.ClusterId)).collect(toList());
            svDataFiltered = true;

            if(!mConfig.ChainIds.isEmpty())
            {
                svData = svData.stream().filter(x -> mConfig.ChainIds.contains(x.ChainId)).collect(toList());
            }
        }

        SvData = svData;

        Fusions = loadFusions(fusionFile).stream().filter(x -> matchOnSampleId(x.SampleId)).collect(toList());
        Exons = VisExons.readExons(geneExonFile).stream().filter(x -> matchOnSampleId(x.SampleId)).collect(toList());
        Segments = VisSegments.readSegments(linksFile).stream().filter(x -> matchOnSampleId(x.SampleId)).collect(toList());

        CopyNumbers = VisCopyNumbers.read(cnaFile).stream().filter(x -> matchOnSampleId(x.SampleId)).collect(toList());

        ProteinDomains = VisProteinDomains.readProteinDomains(proteinFile, Fusions).stream()
                .filter(x -> matchOnSampleId(x.SampleId)).collect(toList());

        if(mConfig.AmberDir != null)
        {
            final String amberBafFile = AmberBAFFile.generateAmberFilenameForReading(mConfig.AmberDir, mConfig.Sample);
            Multimap<Chromosome,AmberBAF> amberBafData = AmberBAFFile.read(amberBafFile, true);
            AmberBAFs.addAll(amberBafData.values());
        }

        if(mConfig.CobaltDir != null)
        {
            final String cobaltRatioFile = CobaltRatioFile.generateFilename(mConfig.CobaltDir, mConfig.Sample);
            ListMultimap<Chromosome, CobaltRatio> cobaltRatiosUnfiltered = CobaltRatioFile.read(cobaltRatioFile);
            List<CobaltRatio> cobaltRatiosFiltered = cobaltRatiosUnfiltered.values().stream().filter(x -> x.tumorGCRatio() != -1).toList();
            CobaltRatios.addAll(cobaltRatiosFiltered);
        }

        if(mConfig.PurpleDir != null)
        {
            final String purpleSegmentFile = PurpleSegment.generateFilename(mConfig.PurpleDir, mConfig.Sample);
            List<PurpleSegment> purpleSegmentsUnfiltered = PurpleSegment.read(purpleSegmentFile);
            List<PurpleSegment> purpleSegmentsFiltered = purpleSegmentsUnfiltered.stream()
                    .filter(x -> x.GermlineState == GermlineStatus.DIPLOID).toList();
            PurpleSegments.addAll(purpleSegmentsFiltered);
        }

        boolean clustersOrChainsProvided = !mConfig.ChainIds.isEmpty() || !mConfig.ClusterIds.isEmpty();
        boolean anyLinxVisDataEmpty = Segments.isEmpty() || SvData.isEmpty() || CopyNumbers.isEmpty();
        boolean anyUpstreamCnvDataExists =  !AmberBAFs.isEmpty() || !CobaltRatios.isEmpty() || !PurpleSegments.isEmpty();

        if(clustersOrChainsProvided && anyLinxVisDataEmpty)
        {
            VIS_LOGGER.error("sample({}) Cannot plot user specified cluster/chain IDs because Linx VIS data was empty", mConfig.Sample);
            System.exit(1);
        }

        if(anyLinxVisDataEmpty)
        {
            if(anyUpstreamCnvDataExists)
            {
                VIS_LOGGER.info("sample({}) Linx VIS data empty, but proceeding to plotting CNV data", mConfig.Sample);
            }
            else
            {
                VIS_LOGGER.info("sample({}) Linx VIS data and (filtered) CNV data empty - no plots to generate", mConfig.Sample);
            }

            return;
        }

        boolean loadSvData = !mConfig.UseCohortFiles
                && (mConfig.PlotClusterGenes || mConfig.RestrictClusterByGene || !mConfig.SpecificRegions.isEmpty());

        if(loadSvData)
        {
            final String svAnnotationsFile = LinxSvAnnotation.generateFilename(mConfig.SampleDataDir, mConfig.Sample, false);
            final String svAnnotationsFileGermline = LinxSvAnnotation.generateFilename(mConfig.SampleDataDir, mConfig.Sample, true);

            List<LinxSvAnnotation> svAnnotations = Files.exists(Paths.get(svAnnotationsFile)) ?
                    LinxSvAnnotation.read(svAnnotationsFile) : LinxSvAnnotation.read(svAnnotationsFileGermline);

            if(svDataFiltered)
            {
                svAnnotations = svAnnotations.stream().filter(x -> SvData.stream().anyMatch(y -> y.SvId == x.svId())).collect(toList());
            }

            if(mConfig.PlotClusterGenes && !mConfig.ClusterIds.isEmpty())
            {
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    if(mConfig.ClusterIds.contains(svAnnotation.clusterId()))
                    {
                        if(!svAnnotation.geneStart().isEmpty())
                            Arrays.stream(svAnnotation.geneStart().split(ITEM_DELIM)).forEach(x -> mConfig.Genes.add(x));

                        if(!svAnnotation.geneEnd().isEmpty())
                            Arrays.stream(svAnnotation.geneEnd().split(ITEM_DELIM)).forEach(x -> mConfig.Genes.add(x));
                    }
                }
            }
            else if(!mConfig.Genes.isEmpty() && mConfig.ClusterIds.isEmpty() && mConfig.RestrictClusterByGene)
            {
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    if(mConfig.Genes.stream().anyMatch(x -> x.equals(svAnnotation.geneStart()) || x.equals(svAnnotation.geneEnd())))
                    {
                        addClusterId(svAnnotation.clusterId());
                    }
                }

                List<LinxDriver> linxDrivers = LinxDriver.read(LinxDriver.generateFilename(mConfig.SampleDataDir, mConfig.Sample));

                for(LinxDriver linxDriver : linxDrivers)
                {
                    if(mConfig.Genes.stream().anyMatch(x -> x.equals(linxDriver.gene())))
                    {
                        addClusterId(linxDriver.clusterId());
                    }
                }
            }

            if(!mConfig.SpecificRegions.isEmpty() && mConfig.ClusterIds.isEmpty())
            {
                // limit to those clusters covered by an SV in the specific regions
                for(LinxSvAnnotation svAnnotation : svAnnotations)
                {
                    addClusterId(svAnnotation.clusterId());
                }
            }
        }
    }

    private boolean matchOnSampleId(final String sampleId)
    {
        return mConfig.UseCohortFiles ? mConfig.Sample.equals(sampleId) : true;
    }

    private void addClusterId(final int clusterId)
    {
        if(!mConfig.ClusterIds.contains(clusterId))
            mConfig.ClusterIds.add(clusterId);
    }

    public Set<Integer> findReportableClusters()
    {
        Set<Integer> clusterIds = Sets.newHashSet();

        Fusions.stream().forEach(x -> clusterIds.add(x.ClusterId));

        if(mConfig.SampleDataDir != null)
        {
            try
            {
                // reportable disruptions
                final List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(mConfig.SampleDataDir, mConfig.Sample));

                final List<Integer> svIds = breakends.stream()
                        .filter(x -> x.reportedDisruption()).map(x -> x.svId()).collect(toList());

                for(Integer svId : svIds)
                {
                    VisSvData svData = SvData.stream().filter(x -> x.SvId == svId).findFirst().orElse(null);
                    if(svData != null)
                         clusterIds.add(svData.ClusterId);
                }

                final List<LinxDriver> drivers = LinxDriver.read(LinxDriver.generateFilename(mConfig.SampleDataDir, mConfig.Sample));
                drivers.stream().filter(x -> x.clusterId() >= 0).forEach(x -> clusterIds.add(x.clusterId()));
            }
            catch(Exception e)
            {
                VIS_LOGGER.error("sample({}) could not read breakends or drivers: {}", mConfig.Sample, e.toString());
            }
        }

        return clusterIds;
    }

    private List<VisFusion> loadFusions(final String fileName) throws IOException
    {
        if(!Files.exists(Paths.get(fileName)))
            return Lists.newArrayList();

        return VisFusion.read(fileName);
    }
}
