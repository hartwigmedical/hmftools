package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.String.format;
import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_CODING_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_MUTATION;
import static com.hartwig.hmftools.common.wisp.CategoryType.OTHER_SV;
import static com.hartwig.hmftools.common.wisp.CategoryType.SUBCLONAL_MUTATION;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_QUALITY_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.NONREPORTABLE_SV_COUNT;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.SUBCLONAL_COUNT;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;
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

    private static final ProbeEvaluator.Criteria PROBE_CRITERIA = new ProbeEvaluator.Criteria(
            SAMPLE_QUALITY_MIN,
            SAMPLE_GC_TARGET, SAMPLE_GC_TOLERANCE);

    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public SampleVariants(final SampleVariantsConfig config, final RefGenomeInterface refGenome, final ProbeGenerator probeGenerator,
            PanelData panelData)
    {
        mConfig = config;
        mRefGenome = refGenome;
        mProbeGenerator = probeGenerator;
        mPanelData = panelData;
    }

    public record ExtraOutput(
            List<Variant> variants
    )
    {
    }

    public ExtraOutput generateProbes()
    {
        LOGGER.info("Generating sample variant probes");

        List<Variant> variants = Lists.newArrayList();

        variants.addAll(loadReferenceVariants());

        checkSampleDirectories(mConfig.purpleDir(), mConfig.linxDir(), mConfig.linxGermlineDir());
        variants.addAll(SomaticMutation.loadSomatics(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(GermlineMutation.loadGermlineMutations(mConfig.sampleId(), mConfig.purpleDir()));
        variants.addAll(StructuralVariant.loadStructuralVariants(mConfig.sampleId(), mConfig.purpleDir(), mConfig.linxDir()));
        variants.addAll(GermlineSv.loadGermlineStructuralVariants(mConfig.sampleId(), mConfig.linxGermlineDir()));

        ProbeGenerationResult result = generateProbes(variants);
        List<Variant> selectedVariants = variants.stream().filter(Variant::isSelected).toList();

        mPanelData.addResult(result);

        LOGGER.info("Done generating sample variant probes");

        return new ExtraOutput(selectedVariants);
    }

    private List<Variant> loadReferenceVariants()
    {
        if(mConfig.referenceVariantsFile() == null)
        {
            return emptyList();
        }
        else
        {
            List<Variant> referenceVariants = ReferenceMutation.loadKnownMutations(mConfig.referenceVariantsFile());
            if(referenceVariants == null)
            {
                throw new RuntimeException("Failed to load reference variants");
            }
            return referenceVariants;
        }
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
        variants.forEach(variant -> variant.setSelectionStatus(SelectionStatus.NOT_SET));
        variants.sort(new VariantComparator());

        ProbeGenerationResult result = new ProbeGenerationResult();

        List<Variant> selectedVariants = Lists.newArrayList();
        ProximateLocations registeredLocations = new ProximateLocations();
        Map<String, Integer> geneDisruptions = Maps.newHashMap();
        int[] typeCounts = new int[CategoryType.values().length];

        result = result.add(generateProbes(variants, selectedVariants, registeredLocations, geneDisruptions, typeCounts, true));
        result = result.add(generateProbes(variants, selectedVariants, registeredLocations, geneDisruptions, typeCounts, false));

        StringJoiner sj = new StringJoiner(" ");
        for(CategoryType type : CategoryType.values())
        {
            sj.add(format("%s=%d", type, typeCounts[type.ordinal()]));
        }
        LOGGER.debug("Selected variant type counts: {}", sj.toString());

        return result;
    }

    private ProbeGenerationResult generateProbes(final List<Variant> variants, List<Variant> selectedVariants,
            ProximateLocations registeredLocations,
            Map<String, Integer> geneDisruptions, int[] typeCounts, boolean firstPass)
    {
        ProbeGenerationResult result = new ProbeGenerationResult();

        for(Variant variant : variants)
        {
            if(selectedVariants.size() >= SAMPLE_PROBES)
            {
                break;
            }

            if(variant.isSelected())
            {
                continue;
            }

            boolean canSelect;

            if(variant.checkFilters())
            {
                if(firstPass)
                {
                    if(!variant.checkAndRegisterGeneLocation(geneDisruptions))
                    {
                        variant.setSelectionStatus(SelectionStatus.GENE_LOCATIONS);
                        canSelect = false;
                    }
                    else if(exceedsMaxByType(variant.categoryType(), OTHER_SV, typeCounts, NONREPORTABLE_SV_COUNT)
                            || exceedsMaxByType(variant.categoryType(), SUBCLONAL_MUTATION, typeCounts, SUBCLONAL_COUNT))
                    {
                        variant.setSelectionStatus(SelectionStatus.EXCEEDS_COUNT);
                        canSelect = false;
                    }
                    else if(!variant.passNonReportableFilters(true))
                    {
                        // This will be checked again on the second pass.
                        variant.setSelectionStatus(SelectionStatus.FILTERED);
                        canSelect = false;
                    }
                    else
                    {
                        canSelect = true;
                    }
                }
                else
                {
                    // If the variant failed the strict filters before, give it a chance to pass the relaxed filters now.
                    // Otherwise, the variant failed other checks and it can't be selected.
                    if(variant.selectionStatus() == SelectionStatus.FILTERED)
                    {
                        canSelect = variant.passNonReportableFilters(false);
                    }
                    else
                    {
                        canSelect = false;
                    }
                }
            }
            else
            {
                canSelect = true;
            }

            if(canSelect)
            {
                if(variant.checkAndRegisterLocation(registeredLocations))
                {
                    result = result.add(generateProbe(variant, selectedVariants, typeCounts));
                }
                else
                {
                    variant.setSelectionStatus(SelectionStatus.PROXIMATE);
                }
            }
        }

        return result;
    }

    private ProbeGenerationResult generateProbe(Variant variant, List<Variant> selectedVariants, int[] typeCounts)
    {
        VariantProbeData probeData = variant.generateProbe(mRefGenome);

        TargetMetadata metadata = createTargetMetadata(variant);
        List<TargetRegion> targetRegions = probeData.regions().stream()
                .map(region -> new TargetRegion(region, metadata))
                .toList();

        // If there's an alt sequence (i.e. not in the ref genome) then always produce a probe to cover it.
        boolean covered = !probeData.hasAltSequence() && probeData.regions().stream().allMatch(mPanelData::isCovered);

        if(covered)
        {
            return new ProbeGenerationResult(emptyList(), targetRegions, emptyList(), emptyList());
        }
        else
        {
            Probe probe;
            if(probeData.sequence() == null)
            {
                probe = mProbeGenerator.mProbeFactory.createProbeFromRegion(requireNonNull(probeData.start()), metadata).orElseThrow();
            }
            else
            {
                probe = mProbeGenerator.mProbeFactory.createProbeFromSequence(probeData.sequence(), metadata).orElseThrow();
            }
            probe = mProbeGenerator.mProbeEvaluator.evaluateProbe(probe, PROBE_CRITERIA);

            if(probe.accepted())
            {
                selectedVariants.add(variant);
                variant.setSelectionStatus(SelectionStatus.SELECTED);
                ++typeCounts[variant.categoryType().ordinal()];
                return new ProbeGenerationResult(List.of(probe), targetRegions, targetRegions, emptyList());
            }
            else
            {
                String rejectionReason = "Probe does not meet criteria " + PROBE_CRITERIA;
                return ProbeGenerationResult.rejectTargets(targetRegions, rejectionReason);
            }
        }
    }

    private static boolean exceedsMaxByType(final CategoryType variantCategory, final CategoryType categoryType, final int[] typeCounts,
            int maxCount)
    {
        if(variantCategory != categoryType)
        {
            return false;
        }
        if(maxCount < 0)
        {
            return true;
        }
        return typeCounts[categoryType.ordinal()] >= maxCount;
    }

    private static class VariantComparator implements Comparator<Variant>
    {
        public int compare(final Variant first, final Variant second)
        {
            // Order by:
            //   1. category
            //   2. reported
            //   3. VCN / copy number

            if(first.categoryType() != second.categoryType())
            {
                return first.categoryType().ordinal() < second.categoryType().ordinal() ? -1 : 1;
            }

            if(first.reported() != second.reported())
            {
                return first.reported() ? -1 : 1;
            }

            if(first.categoryType() == OTHER_MUTATION || first.categoryType() == OTHER_CODING_MUTATION
                    || first.categoryType() == SUBCLONAL_MUTATION)
            {
                // to randomise selection of the lowest priority variants (and avoid multiple selections from highly amplified regions)
                // use the inverse of position as the final comparison
                int locationHash1 = ((SomaticMutation) first).locationHash();
                int locationHash2 = ((SomaticMutation) second).locationHash();

                if(locationHash1 != locationHash2)
                {
                    return locationHash1 < locationHash2 ? -1 : 1;
                }
            }
            else if(first.copyNumber() != second.copyNumber())
            {
                return first.copyNumber() > second.copyNumber() ? -1 : 1;
            }

            return 0;
        }
    }

    private TargetMetadata createTargetMetadata(final Variant variant)
    {
        return new TargetMetadata(TargetMetadata.Type.SAMPLE_VARIANT, variant.description());
    }
}
