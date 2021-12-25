
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - SAGE"

dependencies {
    implementation(project(":hmf-common"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.sage.SageApplication")
}
