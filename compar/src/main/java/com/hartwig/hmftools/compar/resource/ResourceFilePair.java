package com.hartwig.hmftools.compar.resource;

import java.nio.file.Path;

import org.jetbrains.annotations.Nullable;

public record ResourceFilePair(ResourceDefinition definition, String fileName, @Nullable Path oldPath, @Nullable Path newPath)
{
}
