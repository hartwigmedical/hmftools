package com.hartwig.hmftools.finding;

import java.lang.reflect.*;
import java.util.*;

import com.hartwig.hmftools.finding.datamodel.VisualisationFile;

import org.jspecify.annotations.Nullable;

public class VisualisationFileFinder {

    public static List<VisualisationFile> find(@Nullable Object root)
    {
        List<VisualisationFile> result = new ArrayList<>();
        Set<Object> visited = Collections.newSetFromMap(new IdentityHashMap<>());
        visit(root, result, visited);
        return result;
    }

    private static void visit(@Nullable Object obj, List<VisualisationFile> result, Set<Object> visited)
    {
        if (obj == null || visited.contains(obj))
        {
            return;
        }

        visited.add(obj);

        if (obj instanceof VisualisationFile vf)
        {
            result.add(vf);
            return;
        }

        Class<?> clazz = obj.getClass();

        // primitives / simple types
        if (clazz.isPrimitive()
                || clazz.getName().startsWith("java.lang")
                || clazz.isEnum())
        {
            return;
        }

        // collections
        if (obj instanceof Iterable<?> iterable)
        {
            for (Object item : iterable)
            {
                visit(item, result, visited);
            }
            return;
        }

        // maps
        if (obj instanceof Map<?, ?> map) {
            for (Object value : map.values()) {
                visit(value, result, visited);
            }
            return;
        }

        // arrays
        if (clazz.isArray())
        {
            int length = Array.getLength(obj);
            for (int i = 0; i < length; i++)
            {
                visit(Array.get(obj, i), result, visited);
            }
            return;
        }

        // records
        if (clazz.isRecord())
        {
            for (RecordComponent rc : clazz.getRecordComponents())
            {
                try
                {
                    Method accessor = rc.getAccessor();
                    Object value = accessor.invoke(obj);
                    visit(value, result, visited);
                }
                catch (Exception ignored) {}
            }
            return;
        }

        // normal objects
        for (Field field : clazz.getDeclaredFields())
        {
            field.setAccessible(true);
            try
            {
                Object value = field.get(obj);
                visit(value, result, visited);
            }
            catch (Exception ignored) {}
        }
    }
}