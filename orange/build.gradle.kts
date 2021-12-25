
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - ORANGE"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    implementation(libs.log4j.slf4j.impl)
    implementation(libs.kernel)
    implementation(libs.layout)
    implementation(libs.io)
    testImplementation(libs.junit)
    testImplementation(project(path = ":hmf-common", configuration = "tests"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.orange.OrangeApplication")
}
