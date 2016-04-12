#pragma once

#include "dxutil.h"

void SceneInit();

void SceneResize(
    int windowWidth, int windowHeight,
    int renderWidth, int renderHeight);

bool SceneHandleEvent(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

void ScenePaint(ID3D11RenderTargetView* pBackBufferRTV);