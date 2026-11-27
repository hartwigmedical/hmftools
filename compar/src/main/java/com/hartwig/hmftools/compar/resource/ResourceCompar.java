package com.hartwig.hmftools.compar.resource;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;
import static com.hartwig.hmftools.compar.common.MismatchType.DUPLICATE_KEY;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_BOTH;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_NEW;
import static com.hartwig.hmftools.compar.common.MismatchType.INVALID_OLD;
import static com.hartwig.hmftools.compar.common.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.OLD_ONLY;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ResourceCompar
{
    private final ResourceComparConfig mConfig;

    public ResourceCompar(final ConfigBuilder configBuilder)
    {
        mConfig = new ResourceComparConfig(configBuilder);
    }

    public boolean run()
    {
        CMP_LOGGER.info("comparing resources: old({}) new({})", mConfig.OldResourcesDir, mConfig.NewResourcesDir);

        try(ResourceMismatchWriter writer = new ResourceMismatchWriter(mConfig))
        {
            List<Path> oldFiles = findFiles(mConfig.OldResourcesDir);
            List<Path> newFiles = findFiles(mConfig.NewResourcesDir);

            for(ResourceType type : ResourceType.values())
            {
                ResourceDefinition definition = ResourceDefinition.create(type);

                for(ResourceFilePair filePair : pairFiles(definition, oldFiles, newFiles))
                {
                    compareFiles(filePair, writer);
                }
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to compare resources: {}", e.toString());
            return false;
        }

        CMP_LOGGER.info("resource comparison complete");
        return true;
    }

    private static List<Path> findFiles(final String directory) throws IOException
    {
        try(Stream<Path> paths = Files.walk(Paths.get(directory)))
        {
            return paths.filter(Files::isRegularFile).collect(Collectors.toList());
        }
    }

    private static List<ResourceFilePair> pairFiles(
            final ResourceDefinition definition, final List<Path> oldFiles, final List<Path> newFiles) throws IOException
    {
        Map<String,Path> oldFilesByName = filesByName(definition, oldFiles);
        Map<String,Path> newFilesByName = filesByName(definition, newFiles);

        Set<String> fileNames = Sets.newTreeSet(oldFilesByName.keySet());
        fileNames.addAll(newFilesByName.keySet());

        return fileNames.stream()
                .map(x -> new ResourceFilePair(definition, x, oldFilesByName.get(x), newFilesByName.get(x)))
                .collect(Collectors.toList());
    }

    private static Map<String,Path> filesByName(final ResourceDefinition definition, final List<Path> files) throws IOException
    {
        Map<String,Path> filesByName = Maps.newHashMap();

        for(Path path : files)
        {
            String fileName = path.getFileName().toString();

            if(!definition.FileNamePattern.matcher(fileName).matches())
                continue;

            Path otherPath = filesByName.put(fileName, path);

            if(otherPath != null)
                throw new IOException(format("file name(%s) found twice, files cannot be paired: %s and %s", fileName, otherPath, path));
        }

        return filesByName;
    }

    private static void compareFiles(final ResourceFilePair filePair, final ResourceMismatchWriter writer) throws IOException
    {
        ResourceDefinition definition = filePair.definition();
        List<ResourceMismatch> mismatches = Lists.newArrayList();
        ResourceFile oldFile = null;
        ResourceFile newFile = null;

        if(filePair.oldPath() == null || filePair.newPath() == null)
        {
            mismatches.add(new ResourceMismatch(filePair.oldPath() != null ? OLD_ONLY : NEW_ONLY));
        }
        else if(Files.mismatch(filePair.oldPath(), filePair.newPath()) >= 0)
        {
            oldFile = readFile(definition, filePair.oldPath());
            newFile = readFile(definition, filePair.newPath());

            if(oldFile == null && newFile == null)
                mismatches.add(new ResourceMismatch(INVALID_BOTH));
            else if(oldFile == null)
                mismatches.add(new ResourceMismatch(INVALID_OLD));
            else if(newFile == null)
                mismatches.add(new ResourceMismatch(INVALID_NEW));
            else
                mismatches.addAll(compareRows(definition, oldFile, newFile));
        }

        CMP_LOGGER.info("{} file({}) mismatches({})", definition.Type, filePair.fileName(), mismatches.size());
        writer.writeMismatches(filePair, mismatches, oldFile, newFile);
    }

    @Nullable
    private static ResourceFile readFile(final ResourceDefinition definition, final Path path)
    {
        try
        {
            return new ResourceFile(definition, path);
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("{} file({}) is invalid: {}", definition.Type, path, e.getMessage());
            return null;
        }
    }

    private static List<ResourceMismatch> compareRows(
            final ResourceDefinition definition, final ResourceFile oldFile, final ResourceFile newFile)
    {
        List<ResourceMismatch> mismatches = Lists.newArrayList();

        List<String> oldColumns = oldFile.columns();
        List<String> newColumns = newFile.columns();
        String removedColumns = oldColumns.stream().filter(x -> !newColumns.contains(x)).collect(Collectors.joining(ITEM_DELIM));
        String addedColumns = newColumns.stream().filter(x -> !oldColumns.contains(x)).collect(Collectors.joining(ITEM_DELIM));

        if(!removedColumns.isEmpty() || !addedColumns.isEmpty())
        {
            List<String> diffs = List.of(format("columns(%s/%s)", removedColumns, addedColumns));
            mismatches.add(new ResourceMismatch(VALUE, "", diffs, Collections.emptyList(), Collections.emptyList()));
        }

        List<FieldInfo> fields = definition.Fields.stream().map(x -> x.Info).collect(Collectors.toList());
        boolean sameFileColumns = oldFile.fileColumns().equals(newFile.fileColumns());
        Map<String,List<String>> oldLines = oldFile.linesByKey();
        Map<String,List<String>> newLines = newFile.linesByKey();

        for(Map.Entry<String,List<String>> entry : oldLines.entrySet())
        {
            List<String> keyNewLines = newLines.getOrDefault(entry.getKey(), Collections.emptyList());

            ResourceMismatch mismatch = compareKey(
                    fields, sameFileColumns, entry.getKey(), entry.getValue(), oldFile, keyNewLines, newFile);

            if(mismatch != null)
                mismatches.add(mismatch);
        }

        for(Map.Entry<String,List<String>> entry : newLines.entrySet())
        {
            if(!oldLines.containsKey(entry.getKey()))
            {
                mismatches.add(compareKey(
                        fields, sameFileColumns, entry.getKey(), Collections.emptyList(), oldFile, entry.getValue(), newFile));
            }
        }

        return mismatches;
    }

    @Nullable
    private static ResourceMismatch compareKey(
            final List<FieldInfo> fields, final boolean sameFileColumns, final String key, final List<String> oldLines,
            final ResourceFile oldFile, final List<String> newLines, final ResourceFile newFile)
    {
        if(sameFileColumns && oldLines.equals(newLines))
            return null;

        if(oldLines.size() > 1 || newLines.size() > 1)
        {
            List<String> diffs = List.of(format("rows(%d/%d)", oldLines.size(), newLines.size()));
            return new ResourceMismatch(DUPLICATE_KEY, key, diffs, oldLines, newLines);
        }

        if(newLines.isEmpty())
            return new ResourceMismatch(OLD_ONLY, key, Collections.emptyList(), oldLines, newLines);

        if(oldLines.isEmpty())
            return new ResourceMismatch(NEW_ONLY, key, Collections.emptyList(), oldLines, newLines);

        List<String> diffs = findDiffs(fields, oldFile.fieldValues(oldLines.get(0)), newFile.fieldValues(newLines.get(0)));
        return diffs.isEmpty() ? null : new ResourceMismatch(VALUE, key, diffs, oldLines, newLines);
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder("Compar");
        ResourceComparConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        if(!new ResourceCompar(configBuilder).run())
            System.exit(1);
    }
}
