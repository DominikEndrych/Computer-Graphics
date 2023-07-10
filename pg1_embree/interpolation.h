#pragma once
Color3f bilinInterpol(Texture* texture, float x, float y) {
    if (x - 1 < 0 || x + 1 > texture->width() || y - 1 < 0 || y + 1 > texture->height()) {
        return texture->get_texel(x, y);
    }

    Color3f v1 = texture->get_texel(floor(x), floor(y));
    Color3f v2 = texture->get_texel(floor(x), ceil(y));
    Color3f v3 = texture->get_texel(ceil(x), floor(y));
    Color3f v4 = texture->get_texel(ceil(x), ceil(y));

    float tmp;
    float q1 = (1 - modf(x, &tmp)) * (1 - modf(y, &tmp));
    float q2 = (1 - modf(x, &tmp)) * (modf(y, &tmp));
    float q3 = (modf(x, &tmp)) * (1 - modf(y, &tmp));
    float q4 = (modf(x, &tmp)) * (modf(y, &tmp));

    Color3f c1 = Color3f{ v1.r * q1, v1.g * q1, v1.b * q1 };
    Color3f c2 = Color3f{ v2.r * q2, v2.g * q2, v2.b * q2 };
    Color3f c3 = Color3f{ v3.r * q3, v3.g * q3, v3.b * q3 };
    Color3f c4 = Color3f{ v4.r * q4, v4.g * q4, v4.b * q4 };

    Color3f C = Color3f{ c1.r + c2.r + c3.r + c4.r, c1.g + c2.g + c3.g + c4.g, c1.b + c2.b + c3.b + c4.b };

    return C;
}