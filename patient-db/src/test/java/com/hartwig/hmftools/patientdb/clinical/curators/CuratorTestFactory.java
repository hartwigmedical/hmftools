package com.hartwig.hmftools.patientdb.clinical.curators;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidNode;

import org.jetbrains.annotations.NotNull;

public final class CuratorTestFactory
{

    private static final String TUMOR_LOCATION_MAPPING_TSV = Resources.getResource("curators/test_tumor_location_mapping.tsv").getPath();
    private static final String TUMOR_LOCATION_OVERRIDES_TSV =
            Resources.getResource("curators/test_tumor_location_overrides.tsv").getPath();
    private static final String BIOPSY_SITE_MAPPING_TSV = Resources.getResource("curators/test_biopsy_site_mapping.tsv").getPath();
    private static final String TREATMENT_MAPPING_TSV = Resources.getResource("curators/test_treatment_mapping.tsv").getPath();

    private static final String DOID_JSON = Resources.getResource("doid/example_doid.json").getPath();

    private CuratorTestFactory()
    {
    }

    @NotNull
    public static PrimaryTumorCurator primaryTumorCurator()
    {
        try
        {
            List<DoidNode> doidNodes = DiseaseOntology.readDoidOwlEntryFromDoidJson(DOID_JSON).nodes();
            DoidNodesResolver doidNodesResolver = new DoidNodesResolver(doidNodes);
            return new PrimaryTumorCurator(TUMOR_LOCATION_MAPPING_TSV, TUMOR_LOCATION_OVERRIDES_TSV, doidNodesResolver);
        }
        catch(IOException e)
        {
            throw new IllegalStateException("Could not create primary tumor curator!");
        }
    }

    @NotNull
    public static BiopsySiteCurator biopsySiteCurator()
    {
        try
        {
            return new BiopsySiteCurator(BIOPSY_SITE_MAPPING_TSV);
        }
        catch(IOException e)
        {
            throw new IllegalStateException("Could not create biopsy site curator!");
        }
    }

    @NotNull
    public static TreatmentCurator treatmentCurator()
    {
        try
        {
            return new TreatmentCurator(TREATMENT_MAPPING_TSV);
        }
        catch(IOException e)
        {
            throw new IllegalStateException("Could not create treatment curator!");
        }
    }
}
