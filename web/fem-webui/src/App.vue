<script setup lang="ts">
import { RouterLink, RouterView, useRouter, useRoute } from 'vue-router'
import { ElContainer, ElHeader, ElMain, ElButton, ElLink, ElAside, ElMenu } from 'element-plus';
import { Edit } from '@element-plus/icons-vue';

const router = useRouter()
const route = useRoute()

const router_goto = (e: any) => {
  router.push(e)
}
</script>

<template>
  <div class="common-layout">
    <el-container style="height: 100vh;">
      <el-header>
      </el-header>
      <el-container>
        <el-aside width="200px">
          <el-menu @select="router_goto">
            <el-menu-item index="/">
              Main
            </el-menu-item>
            <el-menu-item index="/chart">
              Chart
            </el-menu-item>
          </el-menu>
        </el-aside>
        <el-main>
          <RouterView />
        </el-main>
      </el-container>
    </el-container>
  </div>
</template>

<style scoped>
.common-layout .el-header {
  background-color: var(--el-color-primary-light-7);
  position: relative;
  align-items: center;
}

.common-layout .el-header {
  justify-content: left;
  align-items: center;
}
</style>
