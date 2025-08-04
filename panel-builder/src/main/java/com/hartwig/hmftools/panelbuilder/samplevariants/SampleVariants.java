package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.panelbuilder.PanelData;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SampleVariants
{
    private static final Logger LOGGER = LogManager.getLogger(SampleVariants.class);

    public record ExtraOutput(
            List<Variant> variants
    )
    {
    }

    public static ExtraOutput generateProbes(final SampleVariantsConfig config, final RefGenomeInterface refGenome, PanelData panelData)
    {
        LOGGER.info("Generating sample variant probes");

        List<Variant> commonVariants = new ArrayList<>();

        if(config.referenceVariantsFile() != null)
        {
            List<Variant> referenceVariants = ReferenceMutation.loadKnownMutations(config.referenceVariantsFile());
            if(referenceVariants == null)
            {
                throw new RuntimeException("Failed to load reference variants");
            }
            commonVariants.addAll(referenceVariants);
        }

        List<Variant> variants = Lists.newArrayList();
        variants.addAll(commonVariants);

        checkSampleDirectories(config.purpleDir(), config.linxDir(), config.linxGermlineDir());
        variants.addAll(SomaticMutation.loadSomatics(config.sampleId(), config.purpleDir()));
        variants.addAll(GermlineMutation.loadGermlineMutations(config.sampleId(), config.purpleDir()));
        variants.addAll(StructuralVariant.loadStructuralVariants(config.sampleId(), config.purpleDir(), config.linxDir()));
        variants.addAll(GermlineSv.loadGermlineStructuralVariants(config.sampleId(), config.linxGermlineDir()));

        variants.forEach(variant -> variant.generateSequences(refGenome));

        List<Variant> selectedVariants = VariantSelection.selectVariants(variants);

        // TODO
        //        ProbeGenerationResult result;
        //        panelData.addResult(result);

        LOGGER.info("Done generating sample variant probes");

        return new ExtraOutput(selectedVariants);
    }

    // TODO: needed?
    private static void checkSampleDirectories(final String purpleDir, final String linxDir, final String linxGermlineDir)
    {
        // Allow Linx inputs to be optional.
        if(!Files.exists(Paths.get(purpleDir))
                || (linxDir != null && !Files.exists(Paths.get(linxDir)))
                || (linxGermlineDir != null && !Files.exists(Paths.get(linxGermlineDir))))
        {
            String error = "Missing Purple or Linx directories";
            LOGGER.error(error);
            throw new RuntimeException(error);
        }
    }
}
