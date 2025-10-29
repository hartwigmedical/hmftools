package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.max;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_DRIVER_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAGMENT_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_INDEL_LENGTH_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_NONDRIVER_QUALITY_MIN;
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
import java.util.function.Supplier;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.ProbeEvaluator;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
import com.hartwig.hmftools.panelbuilder.ProbeGenerator;
import com.hartwig.hmftools.panelbuilder.SequenceDefinition;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.UserInputError;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Probes covering variants found in sample data.
// Inputs Linx and Purple data from a previous pipeline run.
// Probes are generated for driver SNV, INDEL, and SV, and then nondriver SNV/INDEL.
// Each variant gets 1 probe which consists of the alt sequence.
public class SampleVariants
{
    private final SampleVariantsConfig mConfig;
    private final ProbeGenerator mProbeGenerator;
    private final PanelData mPanelData;

    private static final ProbeEvaluator.Criteria NONDRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_NONDRIVER_QUALITY_MIN, SAMPLE_NONDRIVER_GC_TARGET, SAMPLE_NONDRIVER_GC_TOLERANCE);
    private static final ProbeEvaluator.Criteria DRIVER_PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_DRIVER_QUALITY_MIN, SAMPLE_DRIVER_GC_TARGET, SAMPLE_DRIVER_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public SampleVariants(final SampleVariantsConfig config, final ProbeGenerator probeGenerator, PanelData panelData)
    {
        mConfig = config;
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
        LOGGER.debug("Selecting up to {} driver variants", max(0, maxProbes));

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
            if(!variant.isDriver())
            {
                continue;
            }

            FilterResult filterResult = driverFilters(variant, geneDisruptions);
            if(filterResult.passed())
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
                registerDisruptedGenes(variant, geneDisruptions);
            }
            else
            {
                LOGGER.trace("Variant failed filter: {{}} reason=\"{}\"", variant, filterResult.failReason());
            }
        }
        return result;
    }

    private ProbeGenerationResult generateNondriverProbes(final List<Variant> variants, int maxProbes)
    {
        LOGGER.debug("Selecting up to {} nondriver variants", max(0, maxProbes));

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

            FilterResult filterResult = nondriverFilters(variant);
            if(filterResult.passed())
            {
                ProbeGenerationResult probeGenResult = generateProbe(variant);
                result = result.add(probeGenResult);
            }
            else
            {
                LOGGER.trace("Variant failed filter: {{}} reason=\"{}\"", variant, filterResult.failReason());
            }
        }
        return result;
    }

    private static FilterResult driverFilters(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        List<Supplier<FilterResult>> filters = new ArrayList<>();

        if(variant instanceof SomaticSv sv)
        {
            if(sv.isReportedDisruption())
            {
                filters.add(() -> vafFilter(sv.vaf()));
                filters.add(() -> tumorFragmentsFilter(sv.tumorFragments()));
            }
        }

        filters.add(() -> geneDisruptionFilter(variant, geneDisruptions));

        return FilterResult.applyFilters(filters);
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

    private static FilterResult nondriverFilters(final SomaticMutation variant)
    {
        return FilterResult.applyFilters(List.of(
                () -> vafFilter(variant.vaf()),
                () -> tumorFragmentsFilter(variant.tumorFragments()),
                () -> indelLengthFilter(variant.indelLength()),
                () -> repeatCountFilter(variant.repeatCount()),
                () -> germlineStatusFilter(variant.germlineStatus())));
    }

    private static FilterResult vafFilter(double vaf)
    {
        return FilterResult.condition(vaf >= SAMPLE_VAF_MIN, "VAF");
    }

    private static FilterResult tumorFragmentsFilter(int tumorFragments)
    {
        return FilterResult.condition(tumorFragments >= SAMPLE_FRAGMENT_COUNT_MIN, "tumor fragments");
    }

    private static FilterResult indelLengthFilter(int indelLength)
    {
        return FilterResult.condition(indelLength <= SAMPLE_INDEL_LENGTH_MAX, "indel length");
    }

    private static FilterResult repeatCountFilter(int repeatCount)
    {
        return FilterResult.condition(repeatCount <= SAMPLE_REPEAT_COUNT_MAX, "repeat count");
    }

    private static FilterResult germlineStatusFilter(final GermlineStatus germlineStatus)
    {
        return FilterResult.condition(germlineStatus == GermlineStatus.DIPLOID, "germline status");
    }

    private static FilterResult geneDisruptionFilter(final Variant variant, final Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            boolean pass = sv.disruptedGenes().stream()
                    .allMatch(gene -> geneDisruptions.getOrDefault(gene, 0) + 1 <= SAMPLE_SV_BREAKENDS_PER_GENE_MAX);
            return FilterResult.condition(pass, "gene disruptions");
        }
        return FilterResult.pass();
    }

    private static void registerDisruptedGenes(final Variant variant, Map<String, Integer> geneDisruptions)
    {
        if(variant instanceof SomaticSv sv)
        {
            sv.disruptedGenes().forEach(gene -> geneDisruptions.put(gene, geneDisruptions.getOrDefault(gene, 0) + 1));
        }
    }

    private ProbeGenerationResult generateProbe(final Variant variant)
    {
        LOGGER.trace("Generating probe for variant: {}", variant);

        SequenceDefinition definition = variant.generateProbe();

        TargetMetadata metadata = createTargetMetadata(variant);
        ProbeEvaluator.Criteria evalCriteria = variant.isDriver() ? DRIVER_PROBE_CRITERIA : NONDRIVER_PROBE_CRITERIA;
        return mProbeGenerator.probe(definition, metadata, evalCriteria, mPanelData);
    }

    private static TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(variant.targetType(), variant.toString());
    }
}
