
plugins {
    id("com.hartwig.java-conventions")
    id("com.hartwig.build-version")
}

description = "HMF Tools - AMBER"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.commons.cli)
    implementation(libs.jcommander)
    implementation(libs.tablesaw.core)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.amber.AmberApplication")
}