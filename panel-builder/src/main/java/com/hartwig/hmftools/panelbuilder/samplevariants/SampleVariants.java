package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INDEL_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NOVEL_SEQUENCE_BASES_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_REPEAT_COUNT_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_SV_BREAKENDS_PER_GENE_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VAF_MIN;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.Probe;
import com.hartwig.hmftools.panelbuilder.ProbeEvaluator;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
import com.hartwig.hmftools.panelbuilder.ProbeGenerator;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.TargetRegion;
import com.hartwig.hmftools.panelbuilder.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleVariants
{
    private final SampleVariantsConfig mConfig;
    private final RefGenomeInterface mRefGenome;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final ProbeEvaluator.Criteria NONDRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_NONDRIVER_QUALITY_MIN, SAMPLE_NONDRIVER_GC_TARGET, SAMPLE_NONDRIVER_GC_TOLERANCE);
    private static final ProbeEvaluator.Criteria DRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_DRIVER_QUALITY_MIN, SAMPLE_DRIVER_GC_TARGET, SAMPLE_DRIVER_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public SampleVariants(final SampleVariantsConfig config, final RefGenomeInterface refGenome, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mProbeGenerator = probeGenerator;
        mPanelData = panelData;
    }

    public void generateProbes()
    {
        LOGGER.info("Generating sample variant probes");

        List<Variant> variants = new ArrayList<>();

        checkSampleDirectories(mConfig.purpleDir(), mConfig.linxDir(), mConfig.linxGermlineDir());
        variants.addAll(SomaticMutation.load(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(GermlineMutation.load(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(SomaticSv.load(mConfig.sampleId(), mConfig.purpleDir(), mConfig.linxDir()));
        variants.addAll(GermlineSv.load(mConfig.sampleId(), mConfig.linxGermlineDir()));

        ProbeGenerationResult result = generateProbes(variants);

        mPanelData.addResult(result);

        LOGGER.info("Done generating sample variant probes");
    }

    private static void checkSampleDirectories(final String purpleDir, final String linxDir, final String linxGermlineDir)
    {
        // Allow Linx inputs to be optional.
        if(!Files.exists(Paths.get(purpleDir))
                || (linxDir != null && !Files.exists(Paths.get(linxDir)))
                || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            throw new UserInputError("Missing Purple or Linx directories");
        }
    }

    public ProbeGenerationResult generateProbes(List<Variant> variants)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        Map<String, Integer> geneDisruptions = new HashMap<>();

        result = result.add(generateDriverProbes(variants, geneDisruptions, mConfig.maxProbes() - result.probes().size()));
        result = result.add(generateNondriverProbes(variants, mConfig.maxProbes() - result.probes().size()));

        return result;
    }

    private ProbeGenerationResult generateDriverProbes(final List<Variant> variants, Map<String, Integer> geneDisruptions, int maxProbes)

    {
        ProbeGenerationResult result = new ProbeGenerationResult();
        for(Variant variant : variants)
        {
            if(result.probes().size() >= maxProbes)
            {
                // We assume that there will be enough probes to cover all the driver variants.
                // If not, then we randomly(?) discarded some drivers, which may be important to handle.
                // Potentially there should be a variant prioritisation scheme to avoid this scenario.
                LOGGER.warn("Filled sample variant probe quota without including all driver variants!");
                // Filled the probe quota.
                break;
            }
            boolean select = variant.isDriver()
                    && driverFilters(variant)
                    && geneDisruptionFilter(variant, geneDisruptions);
            if(select)
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
                registerDisruptedGenes(variant, geneDisruptions);
            }
        }
        return result;
    }

    private static boolean driverFilters(final Variant variant)
    {
        // Only disruption SVs have these filters, other SV types have no additional filters and are always accepted if they are drivers.
        if(variant instanceof SomaticSv sv)
        {
            if(sv.isReportedDisruption())
            {
                return vafFilter(sv.vaf()) && tumorFragmentsFilter(sv.tumorFragments());
            }
        }
        return true;
    }

    private ProbeGenerationResult generateNondriverProbes(final List<Variant> variants, int maxProbes)
    {
        // For nondrivers, we are only interested in somatic SNV/INDEL.
        // Also, some variants are prioritised over others.
        List<SomaticMutation> nondriverVariants = variants.stream()
                .filter(variant -> variant instanceof SomaticMutation)
                .map(variant -> (SomaticMutation) variant)
                .filter(variant -> !variant.isDriver())
                .sorted(new NondriverVariantComparator())
                .toList();

        ProbeGenerationResult result = new ProbeGenerationResult();
        for(SomaticMutation variant : nondriverVariants)
        {
            if(result.probes().size() >= maxProbes)
            {
                // Filled the probe quota.
                break;
            }
            boolean select = !variant.isDriver()
                    && nondriverFilters(variant);
            if(select)
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
            }
        }
        return result;
    }

    private static class NondriverVariantComparator implements Comparator<SomaticMutation>
    {
        @Override
        public int compare(final SomaticMutation v1, final SomaticMutation v2)
        {
            if(v1.isCoding() != v2.isCoding())
            {
                return v1.isCoding() ? -1 : 1;
            }
            if(v1.isClonal() != v2.isClonal())
            {
                return v1.isClonal() ? -1 : 1;
            }
            // Otherwise random but deterministic ordering.
            return Integer.compare(v1.deterministicHash(), v2.deterministicHash());
        }
    }

    private static boolean nondriverFilters(final SomaticMutation variant)
    {
        return vafFilter(variant.vaf())
                && tumorFragmentsFilter(variant.tumorFragments())
                && indelLengthFilter(variant.indelLength())
                && repeatCountFilter(variant.repeatCount())
                && germlineStatusFilter(variant.germlineStatus());
    }

    private static boolean vafFilter(double vaf)
    {
        return vaf >= SAMPLE_VAF_MIN;
    }

    private static boolean tumorFragmentsFilter(int tumorFragments)
    {
        return tumorFragments >= SAMPLE_FRAGMENT_COUNT_MIN;
    }

    private static boolean indelLengthFilter(int indelLength)
    {
        return indelLength <= SAMPLE_INDEL_LENGTH_MAX;
    }

    private static boolean repeatCountFilter(int repeatCount)
    {
        return repeatCount <= SAMPLE_REPEAT_COUNT_MAX;
    }

    private static boolean germlineStatusFilter(final GermlineStatus germlineStatus)
    {
        return germlineStatus == GermlineStatus.DIPLOID;
    }

    private static boolean geneDisruptionFilter(final Variant variant, final Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            return sv.disruptedGenes().stream()
                    .allMatch(gene -> geneDisruptions.getOrDefault(gene, 0) + 1 <= SAMPLE_SV_BREAKENDS_PER_GENE_MAX);
        }
        return true;
    }

    private static void registerDisruptedGenes(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            sv.disruptedGenes().forEach(gene -> geneDisruptions.put(gene, geneDisruptions.getOrDefault(gene, 0)));
        }
    }

    private ProbeGenerationResult generateProbe(final Variant variant)
    {
        LOGGER.trace("Generate probe for variant: {}", variant);

        VariantProbeData probeData = variant.generateProbe(mRefGenome);

        TargetMetadata metadata = createTargetMetadata(variant);
        List<TargetRegion> targetRegions = probeData.regions().stream()
                .map(region -> new TargetRegion(region, metadata))
                .toList();

        // Only do the coverage check for variant where the probe is similar to the ref genome.
        // If the probe is similar (e.g. SNV) then that region could be captured by probing the ref genome sequence,
        // so the variant probe is not needed.
        // If the probe is dissimilar (e.g. large INDEL or SV) then we need the variant probe to capture the variant.
        boolean isNovel = isVariantProbeNovel(probeData);
        boolean covered = !isNovel && probeData.regions().stream().allMatch(mPanelData::isCovered);

        if(covered)
        {
            LOGGER.trace("Variant probe already covered by panel: {}", probeData);
            return ProbeGenerationResult.alreadyCoveredTargets(targetRegions);
        }
        else
        {
            ProbeEvaluator.Criteria evalCriteria = variant.isDriver() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
            Probe probe = mProbeGenerator.mProbeFactory.createProbeFromSequence(probeData.sequence(), metadata).orElseThrow();
            probe = mProbeGenerator.mProbeEvaluator.evaluateProbe(probe, evalCriteria);

            if(probe.accepted())
            {
                return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
            }
            else
            {
                String rejectionReason = "Probe does not meet criteria " + evalCriteria;
                return ProbeGenerationResult.rejectTargets(targetRegions, rejectionReason);
            }
        }
    }

    private static boolean isVariantProbeNovel(final VariantProbeData data)
    {
        ChrBaseRegion start = data.start();
        ChrBaseRegion end = data.end();
        int insertLength = data.insert() == null ? 0 : data.insert().length();
        if(start == null && end == null)
        {
            // Unknown region, assume novel sequence.
            return true;
        }
        else if(start != null && end != null)
        {
            if(start.chromosome().equals(end.chromosome()))
            {
                // SNV, INDEL, or SV on same chromosome.
                if(start.start() > end.start())
                {
                    // Ensure start and end are ordered correctly to calculate the delete length.
                    start = data.end();
                    end = data.start();
                }
                // Clamp to >=0 because theoretically the regions could overlap in the case of an SV.
                int deleteLength = max(end.start() - start.end() - 1, 0);
                int difference = abs(insertLength - deleteLength);
                return difference >= SAMPLE_NOVEL_SEQUENCE_BASES_MIN;
            }
            else
            {
                // SV across different chromosomes.
                return true;
            }
        }
        else
        {
            // Single ended SV.
            return true;
        }
    }

    private static TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(variant.targetType(), variant.toString());
    }
}
