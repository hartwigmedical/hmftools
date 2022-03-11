package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidMetadata {

    private static final Logger LOGGER = LogManager.getLogger(DoidMetadata.class);

    // TODO Consider adding deprecated() and comments() which have been added between 201015 and 220302

    @Nullable
    public abstract DoidDefinition doidDefinition();

    @Nullable
    public abstract List<String> subsets();

    @Nullable
    public abstract List<DoidXref> xrefs();

    @Nullable
    public abstract List<DoidSynonym> synonyms();

    @Nullable
    public abstract List<DoidBasicPropertyValue> basicPropertyValues();

    @Nullable
    @Value.Derived
    public String snomedConceptId() {
        if (xrefs() == null) {
            return null;
        }

        for (DoidXref xref : xrefs()) {
            // Format to look for is SNOMEDCT_US_2020_03_01:109355002
            if (xref.val().contains("SNOMED")) {
                String[] parts = xref.val().split(":");
                if (parts.length == 2 && isLong(parts[1])) {
                    return parts[1];
                } else {
                    LOGGER.warn("Unexpected SNOMED entry found: {}", xref.val());
                }
            }
        }
        return null;
    }

    private static boolean isLong(@NotNull String string) {
        try {
            Long.parseLong(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
}
