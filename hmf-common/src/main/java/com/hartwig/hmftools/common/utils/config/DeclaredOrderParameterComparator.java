package com.hartwig.hmftools.common.utils.config;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterDescription;
import com.beust.jcommander.ParametersDelegate;

import org.jetbrains.annotations.NotNull;

// for use with jcommander, this allows outputting jcommander usage in the order the parameters
// are declared
// Modified from https://groups.google.com/g/jcommander/c/w9vWtzCQfXg written by Brandon Johnson
// Changed to work with ParametersDelegate also
//
public class DeclaredOrderParameterComparator implements Comparator<ParameterDescription> {

    private final List<String> declaredOrderedOptions;

    public DeclaredOrderParameterComparator(@NotNull Class<?> clazz) {
        declaredOrderedOptions = getLongestParameterNamesInDeclaredOrder(clazz);
    }

    @NotNull
    private List<String> getLongestParameterNamesInDeclaredOrder(@NotNull Class<?> clazz) {
        List<String> parameterNames = new ArrayList<>();

        //As of JDK 6.0 fields are returned in the order they were declared.
        Field[] declaredFields = clazz.getDeclaredFields();
        for (Field declaredField : declaredFields) {
            Parameter parameterAnnotation = declaredField.getAnnotation(Parameter.class);
            if (parameterAnnotation != null)
            {
                String[] names = parameterAnnotation.names();
                parameterNames.add(getLongestName(names));
                continue;
            }
            ParametersDelegate delegateAnnotation = declaredField.getDeclaredAnnotation(ParametersDelegate.class);
            if (delegateAnnotation != null)
            {
                parameterNames.addAll(getLongestParameterNamesInDeclaredOrder(declaredField.getType()));
            }
        }

        return parameterNames;
    }

    private String getLongestName(@NotNull String[] names){
        String longestName = "";
        for (String name : names)
        {
            if(name.length() > longestName.length())
                longestName = name;
        }
        return longestName;
    }

    @Override
    public int compare(ParameterDescription o1, ParameterDescription o2) {
        int o1DeclaredIndex = declaredOrderedOptions.indexOf(o1.getLongestName());
        int o2DeclaredIndex = declaredOrderedOptions.indexOf(o2.getLongestName());
        return Integer.compare(o1DeclaredIndex, o2DeclaredIndex);
    }
}
