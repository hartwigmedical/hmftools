
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - SAGE"

dependencies {
    implementation(project(":hmf-common"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.sage.SageApplication")
}
